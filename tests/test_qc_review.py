import json
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from xml.etree.ElementTree import fromstring

from cuttag_web import hospital_pair_review as pair_review
from cuttag_web.hospital_pair import TARGET
from cuttag_web.model import ModelClient, ModelError
from cuttag_web.pipeline import sha256
from cuttag_web.qc_evidence import AUTH_GROUPS, CONTEXT, SCHEMA
from cuttag_web.server import PRESET_ARTIFACTS, Workbench
from cuttag_web.track_review import review_qc_baseline


class LocalModelTransport:
    """Parse the outbound request without contacting any model service."""

    def __init__(self, actions):
        self.actions = iter(actions)
        self.requests = []

    def open(self, request, timeout=None):
        body = json.loads(request.data)
        self.requests.append(json.loads(body['messages'][1]['content']))
        decision = next(self.actions)
        content = decision if isinstance(decision, str) else json.dumps(decision)
        envelope = {'choices': [{'message': {'content': content}}]}
        return FakeResponse(envelope)


class FakeResponse:
    def __init__(self, envelope):
        self.envelope = envelope

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self, limit):
        return json.dumps(self.envelope).encode()


def decision(action, evidence, size=None):
    return {'action': action, 'evidence_used': evidence, 'bin_size': size, 'assessment': '仅用于测试',
            'reason': '测试执行流程', 'missing_information': []}


def pair_fixture():
    return {
        'bin_size': 25, 'peak_settings': {'qvalue': .01, 'mode': 'exploratory_narrow', 'blacklist_filter': 'not_applied'},
        'peak_counts': {'baseline': 80, 'with-igg': 50}, 'fragments': {'target': 1000, 'igg': 500},
        'fragment_overlap': {'target': {'baseline': {'fraction': .2}, 'with-igg': {'fraction': .16}},
                             'igg': {'baseline': {'fraction': .1}, 'with-igg': {'fraction': .01}}},
        'spikein': {'addition_stage': 'after_fragmentation', 'dose_pg': 1, 'cells_per_dose': 200000,
                    'actual_cell_counts': None, 'historical_fragments': {'target': 20, 'igg': 100}},
    }


class QcReviewTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        (self.root / '.local').mkdir()
        (self.root / '.env').write_text('CUTTAG_LLM_BASE_URL=https://model.example/v1\nCUTTAG_LLM_MODEL=fixture\n'
                                       'CUTTAG_LLM_API_KEY=fixture\nCUTTAG_ALLOW_MINIMAL_SUMMARY=1\n')
        core = {'data_source': 'authorized_minimal_summary',
                'parameters': {'region': 'chr1:1000000:2000000', 'normalization': 'CPM', 'bin_size': 25},
                'execution': {'exit_code': 0, 'output_valid': True},
                'bam_counts': {'mapped_segments': 200},
                'bw_statistics': {'interval_count': 3, 'nonzero_intervals': 2}}
        narrow = {'peaks': 7, 'reference_sequences_with_peaks': 1, 'median_width_bp': 125,
                  'summed_width_bp': 850, 'qvalue': .01, 'overlapping_fragments': 30,
                  'total_fragments': 100, 'frip': .3}
        broad = {**narrow, 'peaks': 15, 'median_width_bp': 300, 'broad_cutoff': .1}
        context = {key: sorted(values)[0] for key, values in CONTEXT.items()}
        self.payload = {'data_source': SCHEMA, 'core': core,
                        'qc': {'narrow_peaks': narrow, 'broad_peaks': broad, 'experiment_context': context}}
        evidence = patch('cuttag_web.track_review.make_qc_evidence',
                         return_value=(self.payload, {'narrow_peaks': '/local/private/narrow.bed'}))
        evidence.start()
        self.addCleanup(evidence.stop)

    def authorize(self, base='https://model.example/v1'):
        scope = {'source': 'explicit_user_reply', 'schema': SCHEMA, 'sample_ids': ['fixture-sample'],
                 'base_url': base, 'model': 'fixture', 'allowed': sorted(AUTH_GROUPS)}
        (self.root / '.local/qc_evidence_authorization.json').write_text(json.dumps(scope))

    def review(self, transport):
        model = ModelClient(self.root, opener=transport)
        with patch.object(model, 'status', return_value={'model_connected': True}):
            return review_qc_baseline(self.root, 'resume', 'fixture-sample', None, model, lambda: False, lambda _: None)

    def test_local_preview_and_no_transport_on_unapproved_provider(self):
        transport = LocalModelTransport([])
        self.assertEqual(self.review(transport)['status'], 'evidence_ready')
        self.authorize('https://different.example/v1')
        result = self.review(transport)
        self.assertEqual(result['evidence'], self.payload)
        self.assertEqual(result['status'], 'evidence_ready')
        self.assertEqual(transport.requests, [])

    def test_single_review_igv_reuses_only_its_own_verified_baseline(self):
        workbench = Workbench(self.root)
        baseline_id, review_id = 'run-0123456789ab', 'run-0123456789ac'
        folder = workbench.pipeline.output_dir(baseline_id)
        folder.mkdir(parents=True)
        files = {'coverage.cpm.bw': 'fixture', 'narrow.peaks.bed': 'chr1\t900\t1000\n',
                 'broad.peaks.bed': 'chr1\t1800\t2200\n'}
        artifacts = []
        for name, content in files.items():
            path = folder / name
            path.write_text(content)
            artifacts.append({'name': name, 'bytes': path.stat().st_size})
        workbench.runs[baseline_id] = {'id': baseline_id, 'sample_id': 'fixture-sample', 'mode': 'pipeline',
                                      'status': 'completed_pipeline', 'result': {'validation': 'passed'}, 'artifacts': artifacts}
        review = {'id': review_id, 'sample_id': 'fixture-sample', 'mode': 'review', 'baseline_id': baseline_id,
                  'status': 'completed_review'}
        workbench.runs[review_id] = review
        view, paths = workbench.igv_inputs(review_id)
        self.assertEqual(view['loci'], ['chr1:101-1800', 'chr1:1001-3000'])
        self.assertEqual([track['url'] for track in view['tracks']], [
            f'/api/runs/{baseline_id}/artifacts/{name}' for name in files])
        self.assertEqual(view, workbench.igv_inputs(baseline_id)[0])
        resources = fromstring(workbench.igv_session(review_id)).findall('./Resources/Resource')
        self.assertEqual([node.get('path') for node in resources], [str(path) for path in paths])
        review['sample_id'] = 'different-sample'
        with self.assertRaisesRegex(ValueError, '没有可复用'):
            workbench.igv_inputs(review_id)

    def test_paired_web_review_requires_its_scope_and_reuses_computed_igg(self):
        pair = pair_fixture()
        workbench = Workbench(self.root)
        transport = LocalModelTransport(['模型输出未按 JSON 格式返回', decision('compare_igg', ['baseline']),
                                         decision('keep', ['comparison'])])
        workbench.model = ModelClient(self.root, opener=transport)
        request = {'sample_id': TARGET, 'mode': 'review'}
        self.authorize()
        with self.assertRaisesRegex(ValueError, '外部模型使用范围待确认'):
            workbench.create_run(request)
        self.assertEqual(transport.requests, [])

        scope = {'source': 'explicit_user_reply', 'schema': pair_review.SCHEMA, 'sample_ids': [TARGET],
                 'base_url': 'https://model.example/v1', 'model': 'fixture', 'allowed': sorted(pair_review.AUTH_GROUPS)}
        (self.root / '.local/paired_qc_authorization.json').write_text(json.dumps(scope))
        with patch.object(workbench.hospital_pair, 'sample', return_value={'id': TARGET, 'pair': pair}), \
                patch.object(workbench.model, 'status', return_value={'model_connected': True, 'minimal_summary_allowed': True}):
            with self.assertRaisesRegex(ValueError, '已登记'):
                workbench.create_run({**request, 'baseline_id': 'run-000000000000'})
            run = workbench.create_run(request)
            for _ in range(100):
                run = workbench.get_run(run['id'])
                if run['status'] not in ('queued', 'running'):
                    break
                time.sleep(.01)

        self.assertEqual(run['status'], 'completed_review')
        self.assertEqual(run['actual_model_requests'], 3)
        self.assertEqual(run['track_review']['full_pipeline_executions'], 0)
        self.assertTrue(any('同一份摘要重试' in e['message'] for e in run['events']))
        initial, retry, reassessment = transport.requests
        self.assertEqual(initial, pair_review.make_payload(pair))
        self.assertEqual(retry, initial)
        self.assertEqual(reassessment['comparison']['peak_count'], 50)
        self.assertNotIn(TARGET, json.dumps(transport.requests))
        with self.assertRaisesRegex(ValueError, '字段'):
            workbench.model.review_pair({**initial, 'internal_path': '/data/private'})
        with self.assertRaisesRegex(ValueError, '一次'):
            pair_review.validate_decision(decision('compare_igg', ['comparison']), reassessment)
        invalid = LocalModelTransport(['无 JSON', decision('compare_igg', ['spikein']), decision('keep', ['baseline'])])
        with self.assertRaisesRegex(ModelError, '必须引用'):
            ModelClient(self.root, opener=invalid).review_pair(initial)
        self.assertEqual(invalid.requests, [initial, initial])

    def test_adaptive_review_waits_for_acceptance_and_reassesses_returned_metrics(self):
        pair = pair_fixture()
        group = {'total': 40, 'target_higher': 15, 'igg_higher': 5, 'similar': 20, 'unscored': 0,
                 'target_local_peak': 8, 'median_target_cpm': .3, 'median_igg_cpm': .1}
        pair['igv_screen'] = {'igv_examples_reviewed': 0, 'ratio_threshold': 2,
                              'removed': group, 'overlaps_controlled': group}
        tune = {**decision('tune_igg_qvalue', ['comparison', 'comparison.igv_screen.removed.target_higher']), 'qvalue': .005}
        prefer = {**decision('prefer_rerun', ['comparison']), 'qvalue': None}
        keep = {**decision('keep', ['comparison']), 'qvalue': None}
        compare = {**decision('compare_igg', ['baseline', 'comparison']), 'qvalue': None}
        tuned = {'qvalue': .005, 'peak_count': 38,
                 'target': {'fragments': 1000, 'frip': .14}, 'igg': {'fragments': 500, 'frip': .004}}
        change = {'total': 3, 'target_higher': 2, 'igg_higher': 0, 'similar': 1, 'unscored': 0, 'target_local_peak': 1}
        tuned['igv_changes'] = {'removed_now_overlapping': change, 'overlapping_now_missing': change,
                                'median_igg_cpm_in_rerun_overlap': .08}
        scope = {'source': 'explicit_user_reply', 'schema': pair_review.SCHEMA, 'sample_ids': [TARGET],
                 'base_url': 'https://model.example/v1', 'model': 'fixture', 'allowed': sorted(pair_review.AUTH_GROUPS)}
        (self.root / '.local/paired_qc_authorization.json').write_text(json.dumps(scope))
        workbench = Workbench(self.root)
        transport = LocalModelTransport([tune, prefer, prefer, tune, keep, compare])
        workbench.model = ModelClient(self.root, opener=transport)

        def wait(run_id):
            for _ in range(200):
                run = workbench.get_run(run_id)
                if run['status'] not in ('queued', 'running'):
                    return run
                time.sleep(.01)
            self.fail('任务未结束')

        def synthetic_rerun(run_id, _qvalue, _cancelled, **_options):
            folder = self.root / 'outputs' / 'paired_reruns' / run_id
            folder.mkdir(parents=True)
            (folder / 'tuned.peaks.bed').write_text('chr1\t0\t1\n')
            return tuned, []

        with patch.object(workbench.hospital_pair, 'sample', return_value={'id': TARGET, 'pair': pair}), \
                patch.object(workbench.hospital_pair, 'prepare_screen', return_value=pair['igv_screen']) as prepare, \
                patch.object(workbench.model, 'status', return_value={'model_connected': True, 'minimal_summary_allowed': True}), \
                patch('cuttag_web.server.rerun_peaks', side_effect=synthetic_rerun) as tool, \
                patch('cuttag_web.server.build_comparison_figure') as make_figure:
            request = {'sample_id': TARGET, 'mode': 'paired_agent'}
            with self.assertRaisesRegex(ValueError, '先接受'):
                workbench.create_run({**request, 'defer_rerun': False})
            plan = wait(workbench.create_run(request)['id'])
            self.assertEqual(plan['status'], 'awaiting_rerun')
            self.assertEqual(prepare.call_count, 1)
            events = [event['message'] for event in plan['events']]
            self.assertLess(events.index('本地证据已准备，正在等待 Agent 决定。'), events.index('已发起完整质控摘要请求。'))
            self.assertEqual(plan['track_review']['answer']['decision']['evidence_used'], ['comparison', 'igv_screen'])
            tool.assert_not_called()
            make_figure.assert_not_called()
            resumed = wait(workbench.decide_paired_rerun(plan['id'], True)['id'])
            self.assertEqual(resumed['status'], 'completed_review')
            self.assertTrue(resumed['paired_rerun_accepted'])
            self.assertTrue(any('已接受 Agent' in event['message'] for event in resumed['events']))
            self.assertEqual(resumed['track_review']['reassessment']['decision'], prefer)
            self.assertEqual(resumed['track_review']['peak_executions'], 1)
            self.assertEqual(prepare.call_count, 1)
            tool.assert_called_once()
            self.assertEqual(make_figure.call_args.args[:2], (self.root.resolve(), plan['id']))
            self.assertEqual(tool.call_args.args[:2], (plan['id'], .005))
            self.assertEqual(transport.requests[1]['comparison']['rerun'], tuned)
            with self.assertRaisesRegex(ValueError, '没有可继续'):
                workbench.resume_paired_rerun(plan['id'])
            again = wait(workbench.replay_paired_rerun(resumed['id'])['id'])
            self.assertEqual(again['status'], 'completed_review')
            self.assertNotEqual(again['id'], resumed['id'])
            self.assertEqual(again['reused_plan_id'], resumed['id'])
            self.assertEqual(again['actual_model_requests'], 1)
            self.assertEqual(again['track_review']['reassessment']['decision'], prefer)
            self.assertEqual(again['track_review']['peak_executions'], 1)
            self.assertEqual(tool.call_args.args[:2], (again['id'], .005))
            self.assertEqual(transport.requests[2], transport.requests[1])
            self.assertEqual(workbench.get_run(resumed['id']), resumed)
            pending = wait(workbench.create_run({**request, 'defer_rerun': True})['id'])
            self.assertEqual(pending['status'], 'awaiting_rerun')
            with self.assertRaisesRegex(ValueError, '请选择'):
                workbench.decide_paired_rerun(pending['id'], 'false')
            declined = workbench.decide_paired_rerun(pending['id'], False)
            self.assertEqual(declined['status'], 'completed_review')
            self.assertFalse(declined['paired_rerun_accepted'])
            self.assertEqual(declined['track_review']['peak_executions'], 0)
            self.assertEqual(declined['actual_model_requests'], 1)
            with self.assertRaisesRegex(ValueError, '没有可继续'):
                workbench.decide_paired_rerun(pending['id'], True)
            finished = wait(workbench.create_run(request)['id'])
            self.assertEqual(finished['status'], 'completed_review')
            self.assertEqual(finished['track_review']['peak_executions'], 0)
            with self.assertRaisesRegex(ValueError, '没有已执行'):
                workbench.replay_paired_rerun(finished['id'])
            compared = wait(workbench.create_run(request)['id'])
            self.assertEqual(compared['status'], 'completed_review')
            self.assertTrue(compared['track_review']['existing_comparison_selected'])
            self.assertEqual(compared['track_review']['answer']['decision'], compare)
            self.assertEqual(compared['track_review']['peak_executions'], 0)
            self.assertEqual(compared['actual_model_requests'], 1)
            self.assertEqual(tool.call_count, 2)
            self.assertEqual(make_figure.call_count, 2)
        with self.assertRaisesRegex(ValueError, '允许的 q 值'):
            pair_review.validate_decision({**tune, 'qvalue': .1}, transport.requests[0], adaptive=True)
        with self.assertRaisesRegex(ValueError, '回复格式无效'):
            pair_review.validate_decision(tune, transport.requests[1], adaptive=True)
        self.assertEqual(transport.requests[0]['experiment_context']['igv_review'], 'automated_screened')
        self.assertEqual(transport.requests[0]['comparison']['igv_screen']['removed'], group)
        with self.assertRaisesRegex(ValueError, '已提供'):
            pair_review.validate_decision({**keep, 'evidence_used': ['comparison.igv_screen.unknown']},
                                          transport.requests[0], adaptive=True)
        with self.assertRaisesRegex(ValueError, '待确认范围'):
            pair_review.validate_decision({**keep, 'missing_information': ['igv_review']},
                                          transport.requests[0], adaptive=True)
        with self.assertRaisesRegex(ValueError, '回复格式无效'):
            pair_review.validate_decision({**keep, 'action': 'request_evidence', 'missing_information': ['control_design']},
                                          transport.requests[0], adaptive=True)

    def test_history_archive_preserves_record_and_restores_visibility(self):
        workbench = Workbench(self.root)
        run = {'id': 'run-0123456789ab', 'mode': 'review', 'status': 'completed_review',
               'created_at': '2026-11-26T00:00:00+00:00'}
        workbench.runs[run['id']] = run
        workbench.save(run)
        with self.assertRaises(ValueError):
            workbench.set_run_archived(run['id'], 'false')
        with self.assertRaises(ValueError):
            workbench.set_run_archived('run-000000000000', True)
        self.assertEqual(workbench.set_run_archived(run['id'], True), {'archived_run_ids': [run['id']]})
        self.assertEqual(workbench.get_run(run['id']), run)
        self.assertEqual(Workbench(self.root).archived_run_ids(), [run['id']])
        self.assertEqual(workbench.set_run_archived(run['id'], False), {'archived_run_ids': []})
        self.assertEqual(workbench.get_run(run['id']), run)

    def test_presets_reuse_authorized_payload_and_block_wrong_action(self):
        context = self.payload['qc']['experiment_context']
        context['matched_igg'] = 'unavailable'
        self.payload['qc']['spikein'] = {'spikein_fragments': 100, 'human_fragments': 5000, 'scale_numerator': 10000,
                                        'scale_factor': 100, 'bin_bp': 10, 'same_alignment_processing_as_cpm': False,
                                        'same_bin_as_cpm': False}
        correct = {'action': 'keep_existing_peaks', 'reason': '当前样本没有匹配 IgG'}
        other = {'action': 'generate_pair', 'reason': '有配套 Spike-in 计数'}
        transport = LocalModelTransport([other, correct, other])
        model = ModelClient(self.root, opener=transport)
        with self.assertRaises(ModelError):
            model.review_preset(self.payload, 'fixture-sample', 'igg')
        self.assertEqual(transport.requests, [])
        self.authorize()
        with self.assertRaisesRegex(ModelError, '未执行建议'):
            model.review_preset(self.payload, 'fixture-sample', 'igg')
        self.assertEqual(model.review_preset(self.payload, 'fixture-sample', 'igg')['decision'], correct)
        self.assertEqual(model.review_preset(self.payload, 'fixture-sample', 'normalization')['decision'], other)
        self.assertEqual(transport.requests, [self.payload] * 3)

    def test_web_preset_executes_approved_branch_and_reopens_its_artifacts(self):
        self.authorize()
        (self.root / 'data/private/fixture-sample').mkdir(parents=True)
        workbench = Workbench(self.root)
        review_id, baseline = 'run-0123456789ab', 'run-0123456789ac'
        self.payload['qc']['spikein'] = {'spikein_fragments': 100, 'human_fragments': 5000, 'scale_numerator': 10000,
                                        'scale_factor': 100, 'bin_bp': 10, 'same_alignment_processing_as_cpm': False,
                                        'same_bin_as_cpm': False}
        workbench.save({'id': review_id, 'mode': 'review', 'status': 'completed_review', 'sample_id': 'fixture-sample',
                        'baseline_id': baseline, 'track_review': {'sent_summaries': [self.payload]},
                        'created_at': '2026-11-19T00:00:00+00:00'})
        response = {'action': 'generate_pair', 'reason': '可进行两种标准化的工程对比'}
        transport = LocalModelTransport([response])
        workbench.model = ModelClient(self.root, opener=transport)

        def synthetic_tracks(_, folder, inputs, *, compare, cancelled):
            self.assertTrue(compare)
            self.assertFalse(cancelled())
            folder.mkdir()
            artifacts = []
            for name in PRESET_ARTIFACTS:
                path = folder / name
                path.write_bytes(b'synthetic, not a usable result')
                artifacts.append({'name': name, 'bytes': path.stat().st_size, 'sha256': sha256(path)})
            record = {'status': 'completed_comparison', 'artifacts': artifacts, 'elapsed_seconds': 1,
                      'spikein_fragments': 100, 'scale_factor': 100,
                      'result': {'validation': 'passed', 'cpm': {'bin_size': 25},
                                 'local_spikein_to_cpm_ratio': 0.5, 'compared_preview_bins': 5}}
            (folder / 'record.json').write_text(json.dumps(record))
            return record

        with patch.object(workbench.model, 'status', return_value={'model_connected': True}), \
             patch('tools.run_cuttag_preset.matched_inputs', return_value=(None, None, None, None, 100, 100, 5000)), \
             patch('tools.run_cuttag_preset.run_tracks', side_effect=synthetic_tracks) as runner:
            created = workbench.create_run({'mode': 'preset', 'preset': 'normalization',
                                            'sample_id': 'fixture-sample', 'review_id': review_id})
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline:
                result = workbench.get_run(created['id'])
                if result['status'] not in ('queued', 'running'):
                    break
                time.sleep(.01)
            runner.assert_called_once()
        self.assertEqual(result['status'], 'completed_comparison')
        self.assertEqual(result['answer']['decision'], response)
        self.assertEqual(result['coverage_executions'], 2)
        self.assertEqual(transport.requests, [self.payload])
        reopened = Workbench(self.root)
        self.assertEqual(reopened.get_run(created['id'])['source'], 'web')
        self.assertEqual(reopened.artifact(created['id'], 'coverage.cpm.bw')[0].name, 'coverage.cpm.bw')
        with self.assertRaises(ValueError):
            reopened.artifact(created['id'], 'analysis.bam')

    def test_model_chooses_peak_comparison_and_reassesses_its_actual_output(self):
        self.authorize()
        transport = LocalModelTransport([decision('compare_peak_sets', ['narrow_peaks', 'broad_peaks']),
                                         decision('keep', ['comparison', 'narrow_peaks'])])
        result = self.review(transport)
        self.assertEqual(result['status'], 'completed_review')
        self.assertEqual(result['model_requests_attempted'], 2)
        self.assertEqual(result['comparison']['broad_minus_narrow_peaks'], 8)
        self.assertEqual(transport.requests, [self.payload, {**self.payload, 'comparison': result['comparison']}])
        self.assertEqual(result['sent_summaries'], transport.requests)
        self.assertNotIn('/local/private', json.dumps(transport.requests))

    def test_invalid_action_does_not_trigger_comparison(self):
        self.authorize()
        transport = LocalModelTransport([decision('compare_peak_sets', ['narrow_peaks'])])
        with patch('cuttag_web.track_review.compare_evidence') as compare:
            result = self.review(transport)
        self.assertEqual(result['status'], 'failed')
        self.assertNotIn('comparison', result)
        compare.assert_not_called()

    def test_only_unconfirmed_background_can_be_requested(self):
        self.authorize()
        context = self.payload['qc']['experiment_context']
        context['matched_igg'] = 'unavailable'
        request = decision('request_evidence', ['experiment_context'])
        for fields in (['matched_igg'], ['请提供匹配 IgG'], ['target_expectation', 'target_expectation']):
            with self.subTest(fields=fields):
                result = self.review(LocalModelTransport([{**request, 'missing_information': fields}]))
                self.assertEqual(result['status'], 'failed')
                self.assertNotIn('answer', result)
        accepted = {**request, 'missing_information': ['target_expectation']}
        result = self.review(LocalModelTransport([accepted]))
        self.assertEqual(result['status'], 'completed_review')
        self.assertEqual(result['answer']['decision'], accepted)
        context['matched_igg'] = 'unknown'
        accepted = {**request, 'missing_information': ['matched_igg']}
        result = self.review(LocalModelTransport([accepted]))
        self.assertEqual(result['status'], 'completed_review')
        self.assertEqual(result['answer']['decision'], accepted)

    def test_private_metadata_cannot_be_added_to_export(self):
        self.authorize()
        payload = {**self.payload, 'qc': {**self.payload['qc'], 'sample_id': 'fixture-sample'}}
        transport = LocalModelTransport([])
        with self.assertRaises(ValueError):
            ModelClient(self.root, opener=transport).review_qc(payload, 'fixture-sample')
        with self.assertRaises(ModelError):
            ModelClient(self.root, opener=transport).review_qc(self.payload, 'other-sample')
        self.assertEqual(transport.requests, [])
