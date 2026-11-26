"""Synthetic control/API tests; never use private sequences or call a model."""
import http.client
import json
import tempfile
import threading
import time
import unittest
from http.server import ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from cuttag_web.pipeline import DISK_RESERVE, RESUME_START_FREE, PipelineRunner
from cuttag_web.server import Workbench, PROJECT_ROOT, handler_for


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.sample = self.root / 'data/private/synthetic'
        self.support = self.root / 'data/support/synthetic'
        for folder in ('processed', 'raw', 'qc', 'provenance'):
            (self.sample / folder).mkdir(parents=True)
        for name in ('synthetic.hg38.dedup.bam', 'synthetic.hg38.dedup.bam.bai'):
            (self.sample / 'processed' / name).write_bytes(b'synthetic fixture, not a real BAM')
        for i in (1, 2):
            (self.sample / 'raw' / ('synthetic_%d.fq.gz' % i)).write_bytes(b'synthetic fixture')
        for name in ('reference/hg38/hg38.bam-compatible.chrom.sizes', 'annotation/gencode.v44.genes.bam-compatible.bed',
                     'provenance/reference_compatibility.json'):
            path = self.support / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({'status': 'passed'}) if name.endswith('.json') else 'synthetic\n')
        index = self.support / 'reference/hg38/bowtie2'
        index.mkdir()
        for suffix in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
            (index / ('hg38_noalt.' + suffix + '.bt2')).write_bytes(b'synthetic')
        self.runner = PipelineRunner(self.root)
        self.workbench = Workbench(self.root, pipeline_runner=self.runner)

    def tearDown(self):
        self.tmp.cleanup()

    def wait_done(self, run_id):
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            run = self.workbench.get_run(run_id)
            if run['status'] not in ('queued', 'running'):
                return run
            time.sleep(.01)
        self.fail('Synthetic control did not finish')

    def test_readonly_mounts_and_source_whitelist(self):
        sample, support, inputs = self.runner.inputs('synthetic', 'fastq')
        self.assertEqual(len(inputs), 10)
        argv = self.runner.command('run-123456789abc', sample, support, self.root / 'output', 'sha256:synthetic', 'fastq')
        self.assertIn('type=bind,src=' + str(sample) + ',dst=/input,readonly', argv)
        self.assertIn('type=bind,src=' + str(support) + ',dst=/support,readonly', argv)
        self.assertEqual(argv[argv.index('--network') + 1], 'none')
        self.assertEqual(argv[argv.index('--memory') + 1], '6g')
        self.assertNotIn('.env', ' '.join(argv))
        for source in ('../file', 'spikein', 'sh', None):
            with self.assertRaises(ValueError):
                self.runner.inputs('synthetic', source)
        for extra in ({'command': 'sh'}, {'source': 'arbitrary'}, {'target_bin_size': 50}):
            with self.assertRaises(ValueError):
                self.workbench.create_run({'mode': 'pipeline', 'sample_id': 'synthetic', **extra})
        with self.assertRaises(ValueError):
            self.workbench.create_run({'mode': 'demo', 'sample_id': 'synthetic', 'source': 'bam'})

    def test_outside_reference_symlinks_rejected(self):
        path = self.support / 'annotation/gencode.v44.genes.bam-compatible.bed'
        path.unlink()
        outside = self.root / 'outside.bed'
        outside.write_text('synthetic')
        path.symlink_to(outside)
        with self.assertRaises(ValueError):
            self.runner.inputs('synthetic', 'bam')
        for name in ('../escape', 'run-123456789abc/other'):
            with self.assertRaises(ValueError):
                self.runner.output_dir(name)

    def test_raw_qc_metrics_do_not_claim_filtered_bam_mapping_rate(self):
        folder = self.root / 'metrics'
        folder.mkdir()
        (folder / 'fastp.json').write_text(json.dumps({'summary': {'before_filtering': {'total_reads': 100, 'q30_rate': .95},
                                                                                'after_filtering': {'total_reads': 96}}}))
        (folder / 'alignment.flagstat.json').write_text(json.dumps({'QC-passed reads': {'primary': 96, 'primary mapped': 90},
                                                                   'QC-failed reads': {'primary': 0, 'primary mapped': 0}}))
        metrics = self.runner.report_metrics(folder, 'fastq')
        self.assertEqual(metrics['primary_mapping_rate'], 90 / 96)
        self.assertEqual(metrics['raw_reads'], 100)
        self.assertEqual(self.runner.report_metrics(folder, 'bam'), {})

    def test_incomplete_reference_and_output_rejected(self):
        path = self.support / 'provenance/reference_compatibility.json'
        path.write_text(json.dumps({'status': 'failed'}))
        with self.assertRaises(ValueError):
            self.runner.inputs('synthetic', 'bam')
        with self.assertRaises(ValueError):
            self.runner.validate(self.root, 'fastq')

    def test_timeout_stops_container_without_publishing_results(self):
        self.runner.timeout = 0
        process = Mock(returncode=None)
        process.poll.side_effect = lambda: process.returncode
        def stop(_):
            process.returncode = 137
        with patch.object(self.runner, 'probe', return_value={'ready': True, 'image_id': 'synthetic'}), \
             patch('cuttag_web.pipeline.subprocess.Popen', return_value=process), \
             patch.object(self.runner, 'cancel', side_effect=stop):
            result = self.runner.run('run-123456789abc', 'synthetic', 'bam', lambda: False, lambda *_: None)
        self.assertEqual(result['status'], 'timed_out')
        self.assertNotIn('result', result)
        self.assertTrue(result['execution']['inputs_unchanged'])
        self.assertFalse(any(a['name'].endswith('.bw') for a in result['artifacts']))

    def test_resume_disk_guard_preserves_checkpoint_and_cleans_scratch(self):
        previous_id = 'run-0123456789ab'
        folder = self.runner.output_dir(previous_id)
        folder.mkdir(parents=True)
        pair = [folder / name for name in ('analysis.bam', 'analysis.bam.bai')]
        for path in pair:
            path.write_bytes(b'synthetic checkpoint, not a real BAM')
        record = self.root / 'outputs' / 'web_runs' / (previous_id + '.json')
        record.write_text(json.dumps(dict(id=previous_id, mode='pipeline', source='fastq', status='failed', sample_id='synthetic')))
        audit = folder / 'commands.json'
        old_check = [{'argv': ['samtools', name, '/output/analysis.bam'], 'exit_code': 0} for name in ('index', 'quickcheck', 'idxstats')]
        audit.write_text(json.dumps(old_check[:-1]))
        with self.assertRaises(ValueError):
            self.runner.resume_inputs(previous_id)
        audit.write_text(json.dumps(old_check))

        process = Mock(returncode=None)
        process.poll.side_effect = lambda: process.returncode
        run_id = 'run-123456789abc'
        new_folder = self.runner.output_dir(run_id)

        def start(*_, **__):
            (new_folder / '.scratch' / 'deepTools.tmp').write_bytes(b'temporary file')
            return process

        def stop(_):
            process.returncode = 143

        with patch.object(self.runner, 'probe', return_value={'ready': True, 'image_id': 'synthetic'}), \
             patch('cuttag_web.pipeline.shutil.disk_usage', side_effect=[SimpleNamespace(free=RESUME_START_FREE), SimpleNamespace(free=DISK_RESERVE)]), \
             patch('cuttag_web.pipeline.subprocess.Popen', side_effect=start) as launch, \
             patch.object(self.runner, 'cancel', side_effect=stop):
            result = self.runner.run(run_id, 'synthetic', 'resume', lambda: False, lambda *_: None, resume_from=previous_id)
        self.assertEqual(result['status'], 'failed')
        self.assertIn('50 GB', result['execution']['error'])
        self.assertNotIn('result', result)
        argv = launch.call_args.args[0]
        self.assertIn('type=bind,src=' + str(pair[0]) + ',dst=/resume/analysis.bam,readonly', argv)
        self.assertIn('TMPDIR=/output/.scratch', argv)
        self.assertFalse(any('dst=/input' in arg for arg in argv))
        self.assertTrue(result['execution']['inputs_unchanged'])
        self.assertFalse((new_folder / '.scratch').exists())

    def test_single_real_job_across_modes_and_cancellation(self):
        entered = threading.Event()
        def wait(run_id, sample, source, cancelled, phase):
            phase('coverage', 'synthetic control', {})
            entered.set()
            while not cancelled():
                time.sleep(.01)
            return {'status': 'cancelled', 'actual_model_requests': 0}
        with patch.object(self.runner, 'probe', return_value={'ready': True}), \
             patch.object(self.runner, 'run', side_effect=wait), \
             patch.object(self.workbench.coverage, 'probe', return_value={'ready': True}):
            run = self.workbench.create_run({'mode': 'pipeline', 'source': 'bam', 'sample_id': 'synthetic'})
            self.assertTrue(entered.wait(1))
            for mode in ('pipeline', 'coverage'):
                with self.assertRaises(ValueError):
                    self.workbench.create_run({'mode': mode, 'sample_id': 'synthetic'})
            self.workbench.cancel_run(run['id'])
            final = self.wait_done(run['id'])
        self.assertEqual(final['status'], 'cancelled')
        self.assertNotIn('result', final)
        self.assertEqual(final['actual_model_requests'], 0)

    def test_tool_failure_does_not_publish_output(self):
        with patch.object(self.runner, 'probe', return_value={'ready': True}), \
             patch.object(self.runner, 'run', side_effect=ValueError('Synthetic tool failure')):
            run = self.workbench.create_run({'mode': 'pipeline', 'sample_id': 'synthetic'})
            final = self.wait_done(run['id'])
        self.assertEqual(final['status'], 'failed')
        self.assertIn('Synthetic', final['error'])
        with self.assertRaises(ValueError):
            self.workbench.artifact(run['id'], 'coverage.cpm.bw')

    def test_review_reuses_only_matching_baseline_and_saves_model_decision(self):
        baseline = 'run-0123456789ab'
        self.workbench.runs[baseline] = {'id': baseline, 'sample_id': 'synthetic', 'mode': 'pipeline',
                                         'source': 'resume', 'status': 'completed_pipeline'}
        payload = {'mode': 'review', 'sample_id': 'synthetic', 'baseline_id': baseline}
        reviewed = {'model_requests_attempted': 1, 'status': 'completed_review', 'answer': {'decision': {'action': 'keep'}}}
        with patch.object(self.workbench, 'review_sample_ids', return_value=['synthetic']), \
             patch.object(self.workbench.model, 'status', return_value={'model_connected': True, 'minimal_summary_allowed': True}), \
             patch.object(self.runner, 'run') as pipeline, patch('cuttag_web.server.review_baseline', return_value=reviewed) as review:
            with self.assertRaises(ValueError):
                self.workbench.create_run({**payload, 'sample_id': 'other'})
            with self.assertRaises(ValueError):
                self.workbench.create_run({**payload, 'target_bin_size': 50})
            run = self.workbench.create_run(payload)
            final = self.wait_done(run['id'])
        pipeline.assert_not_called()
        self.assertEqual(review.call_args.args[:2], (self.runner.output_dir(baseline), 'resume'))
        self.assertEqual(final['track_review']['answer']['decision']['action'], 'keep')
        self.assertEqual(final['status'], 'completed_review')
        self.assertEqual(final['actual_model_requests'], 1)

    def test_full_run_keeps_baseline_when_optional_model_review_fails(self):
        base = {'status': 'completed_pipeline', 'result': {'validation': 'passed'}, 'artifacts': [{'name': 'coverage.cpm.bw'}]}
        with patch.object(self.workbench, 'review_sample_ids', return_value=['synthetic']), \
             patch.object(self.workbench.model, 'status', return_value={'model_connected': True, 'minimal_summary_allowed': True}), \
             patch.object(self.runner, 'probe', return_value={'ready': True}), patch.object(self.runner, 'run', return_value=base), \
             patch('cuttag_web.server.review_baseline', return_value={'status': 'failed', 'error': '模拟回复错误', 'model_requests_attempted': 1}) as review:
            run = self.workbench.create_run({'mode': 'pipeline', 'sample_id': 'synthetic', 'source': 'fastq', 'model_review': True})
            final = self.wait_done(run['id'])
        self.assertEqual(review.call_args.args[:2], (self.runner.output_dir(run['id']), 'fastq'))
        self.assertEqual(final['status'], 'completed_pipeline')
        self.assertEqual(final['artifacts'], base['artifacts'])
        self.assertEqual(final['track_review']['status'], 'failed')

    def test_artifact_stream_range_head_and_access_limits(self):
        run_id = 'run-123456789abc'
        folder = self.runner.output_dir(run_id)
        folder.mkdir(parents=True)
        content = b'synthetic range test content, not an analysis result'
        (folder / 'tool.log').write_bytes(content)
        self.workbench.runs[run_id] = {'id': run_id, 'mode': 'pipeline', 'status': 'failed',
            'artifacts': [{'name': 'tool.log', 'bytes': len(content)}]}
        server = ThreadingHTTPServer(('127.0.0.1', 0), handler_for(self.workbench, PROJECT_ROOT / 'web/dist'))
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        def request(name='tool.log', method='GET', range_value=None):
            c = http.client.HTTPConnection('127.0.0.1', server.server_port)
            c.request(method, '/api/runs/' + run_id + '/artifacts/' + name, headers={'Range': range_value} if range_value else {})
            response = c.getresponse()
            result = response.status, dict(response.getheaders()), response.read()
            c.close()
            return result
        try:
            self.assertEqual(request()[2], content)
            status, headers, body = request(range_value='bytes=2-7')
            self.assertEqual((status, body), (206, content[2:8]))
            self.assertEqual(headers['Content-Range'], 'bytes 2-7/' + str(len(content)))
            self.assertEqual(request(range_value='bytes=-4')[2], content[-4:])
            self.assertEqual(request(range_value='bytes=500-')[0], 416)
            self.assertEqual(request(method='HEAD')[2], b'')
            self.assertEqual(request(name='coverage.cpm.bw')[0], 404)
            self.assertEqual(request(name='.env')[0], 404)
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


    def test_stopped_pipeline_keeps_written_logs_without_exposing_partial_track(self):
        run_id = 'run-123456789abc'
        folder = self.runner.output_dir(run_id)
        folder.mkdir(parents=True)
        log = folder / 'tool.log'
        log.write_text('synthetic tool started', encoding='utf-8')
        (folder / 'coverage.cpm.bw').write_bytes(b'synthetic partial file')
        self.workbench.runs[run_id] = {'id': run_id, 'mode': 'pipeline', 'status': 'cancelled', 'artifacts': []}
        result = self.workbench.get_run(run_id)
        self.assertEqual([x['name'] for x in result['artifacts']], ['tool.log'])
        self.assertEqual(self.workbench.artifact(run_id, 'tool.log')[0], log)
        with self.assertRaises(ValueError):
            self.workbench.artifact(run_id, 'coverage.cpm.bw')

    def test_older_completed_analysis_survives_more_than_a_hundred_new_reviews(self):
        baseline = 'run-000000000000'
        self.workbench.save({'id': baseline, 'mode': 'pipeline', 'sample_id': 'synthetic', 'status': 'completed_pipeline',
                             'report_metrics': {}, 'created_at': '2026-11-17T00:00:00+00:00'})
        for i in range(1, 102):
            self.workbench.save({'id': f'run-{i:012x}', 'mode': 'review', 'sample_id': 'synthetic', 'status': 'completed_review',
                                 'baseline_id': baseline, 'created_at': '2026-11-18T00:00:00+00:00'})
        reopened = Workbench(self.root, pipeline_runner=self.runner)
        self.assertEqual(reopened.get_run(baseline)['status'], 'completed_pipeline')
        self.assertIn(baseline, [r['id'] for r in reopened.list_runs()])


if __name__ == '__main__':
    unittest.main()
