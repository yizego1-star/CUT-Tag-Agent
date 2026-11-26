"""Protocol and controller checks use synthetic values and do not call an external API."""
import copy
import io
import json
import tempfile
import unittest
import urllib.error
from pathlib import Path

from cuttag_web.agent import drive_agent
from cuttag_web.model import ModelClient, ModelError, minimal_summary, validate_decision, validate_track_review


def decision(action="regenerate_coverage", size=50):
    return {"action": action, "bin_size": size, "assessment": "合成工程检查", "reason": "合成测试目标", "missing_information": []}


def output(size):
    return {"status": "completed", "actual_tool_executions": 0,
            "result": {"bin_size": size, "preview": {"points": [1, 2]}, "validation": "passed", "interval_count": 2, "nonzero_intervals": 2,
                       "sample_id": "SHOULD_NOT_EXPORT", "path": "/private/SHOULD_NOT_EXPORT"},
            "bam_stats": {"mapped_segments": 10, "unmapped_segments": 0, "reference_count": 1, "path": "SHOULD_NOT_EXPORT"}}


class FakeCoverage:
    def __init__(self): self.calls = []
    def run(self, *args, bin_size, round_index):
        self.calls.append((bin_size, round_index))
        return output(bin_size)


class FakeModel:
    def __init__(self, responses): self.responses = iter(responses); self.calls = []
    def decide(self, summary, on_request=None):
        self.calls.append(copy.deepcopy(summary))
        if on_request: on_request()
        value = next(self.responses)
        if isinstance(value, Exception): raise value
        return {"model": "synthetic", "decision": value, "usage": {}}


class ProtocolTest(unittest.TestCase):
    def test_invalid_actions_parameters_and_limit(self):
        for value in (decision("shell"), decision(size=True), decision(size=25), decision(size=100), {**decision(), "command": "forbidden"}, decision("finish", 50)):
            with self.subTest(value=value), self.assertRaises(ModelError): validate_decision(value, 25, 1)
        with self.assertRaises(ModelError): validate_decision(decision(), 25, 0)
        self.assertEqual(validate_decision(decision("finish", None), 50, 0)["action"], "finish")

    def test_summary_exports_only_allowed_fields(self):
        value = output(25)
        summary = minimal_summary(value["result"], value["bam_stats"], 0)
        self.assertNotIn("SHOULD_NOT_EXPORT", json.dumps(summary))
        self.assertEqual(summary["bam_counts"]["mapped_segments"], 10)

    def test_http_contract_private_key_and_error_redaction(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / ".env").write_text("CUTTAG_LLM_BASE_URL=https://example.invalid/v1\nCUTTAG_LLM_API_KEY=synthetic-secret\nCUTTAG_LLM_MODEL=synthetic\nCUTTAG_ALLOW_MINIMAL_SUMMARY=1\n")
            requests = []
            class Opener:
                def open(self, request, **kwargs):
                    requests.append(request)
                    return io.BytesIO(json.dumps({"choices": [{"message": {"content": json.dumps(decision())}}], "usage": {"total_tokens": 10}}).encode())
            client = ModelClient(root, opener=Opener())
            receipt = client.test()
            self.assertTrue(receipt["ok"])
            self.assertNotIn("synthetic-secret", json.dumps(client.status()))
            self.assertNotIn("synthetic-secret", client.receipt_path.read_text())
            self.assertEqual(requests[0].get_header("Authorization"), "Bearer synthetic-secret")
            self.assertNotIn("synthetic-secret", requests[0].data.decode())
            summary = minimal_summary(output(25)["result"], output(25)["bam_stats"], 0)
            for bad in ({**summary, "sample_id": "secret"}, {**summary, "bam_counts": {"path": "secret"}}, {**summary, "round_index": "secret"}):
                with self.assertRaises(ModelError): client.decide(bad)
            self.assertEqual(len(requests), 1)
            class Failure:
                def open(self, *args, **kwargs):
                    raise urllib.error.HTTPError("https://example.invalid", 401, "synthetic-secret", {}, io.BytesIO(b"synthetic-secret"))
            failed = ModelClient(root, opener=Failure()).test()
            self.assertFalse(failed["ok"])
            self.assertNotIn("synthetic-secret", json.dumps(failed))
            self.assertIn("401", failed["error"])

    def test_real_summary_blocked_without_permission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / ".env").write_text("CUTTAG_LLM_BASE_URL=https://example.invalid/v1\nCUTTAG_LLM_API_KEY=synthetic-secret\nCUTTAG_LLM_MODEL=synthetic\n")
            summary = minimal_summary(output(25)["result"], output(25)["bam_stats"], 0)
            with self.assertRaisesRegex(ModelError, "范围尚未确认"): ModelClient(root).decide(summary)

    def test_track_review_has_no_forced_target_and_can_decline_recalculation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / '.env').write_text('CUTTAG_LLM_BASE_URL=https://example.invalid/v1\nCUTTAG_LLM_API_KEY=synthetic-secret\nCUTTAG_LLM_MODEL=synthetic\n')
            accepted = {**decision('request_evidence', None), 'missing_information': ['合成测试：对照信息']}
            payloads = []

            class Opener:
                def open(self, request, **kwargs):
                    payloads.append(json.loads(request.data))
                    final_text = '<think>合成测试包装文本</think>\n```JSON\n' + json.dumps(accepted) + '\n```'
                    return io.BytesIO(json.dumps({'choices': [{'message': {'content': final_text}}]}).encode())

            client = ModelClient(root, opener=Opener())
            summary = {'data_source': 'synthetic', 'parameters': {'region': 'chr1:1000000:2000000', 'bin_size': 25, 'normalization': 'CPM'},
                       'execution': {'exit_code': 0, 'output_valid': True},
                       'bam_counts': {'pre_filter_primary_segments': 100, 'pre_filter_primary_mapped_segments': 95},
                       'bw_statistics': {'preview_mean': 1.5, 'preview_median': 0.5, 'preview_p95': 3.0}}
            answer = client.review_track(summary)
            self.assertEqual(answer['decision']['action'], 'request_evidence')
            self.assertNotIn('requested_bin_size', json.dumps(payloads))
            sent = json.loads(payloads[0]['messages'][1]['content'])
            self.assertEqual(sent['bam_counts']['pre_filter_primary_segments'], 100)
            with self.assertRaises(ModelError):
                client.review_track({**summary, 'goal': {'requested_bin_size': 50}})
            with self.assertRaises(ModelError):
                client.review_track({**summary, 'bw_statistics': {'preview_mean': 1.5, 'path': '/private/file'}})
            with self.assertRaises(ModelError):
                client.review_track({**summary, 'bam_counts': {'frip': 1}})
            self.assertEqual(len(payloads), 1)
        for candidate in (decision('keep', None), decision('compare_resolution', 50), decision('compare_resolution', 100)):
            self.assertEqual(validate_track_review(candidate, 25), candidate)
        for bad in (decision('request_evidence', None), decision('compare_resolution', 10), decision('shell', None)):
            with self.assertRaises(ModelError):
                validate_track_review(bad, 25)

    def test_agent_reexecutes_once_and_reviews(self):
        coverage, model = FakeCoverage(), FakeModel([decision(), decision("finish", None)])
        result = drive_agent("synthetic-run", "synthetic", coverage, model, lambda: False, lambda *_: None)
        self.assertEqual(result["status"], "completed_agent")
        self.assertEqual(coverage.calls, [(25, 0), (50, 1)])
        self.assertEqual(len(model.calls), 2)
        self.assertEqual(model.calls[-1]["reruns_left"], 0)
        self.assertEqual(result["reruns_completed"], 1)
        self.assertEqual(result["actual_tool_executions"], 0, "Synthetic executor must not count as real tool execution")

    def test_bad_advice_and_model_failure_do_not_trigger_rerun(self):
        for response in (decision("shell"), ModelError("合成服务超时")):
            coverage, model = FakeCoverage(), FakeModel([response])
            result = drive_agent("synthetic-run", "synthetic", coverage, model, lambda: False, lambda *_: None)
            self.assertEqual(result["status"], "failed")
            self.assertEqual(len(coverage.calls), 1)
        coverage, model = FakeCoverage(), FakeModel([decision(), decision(size=25)])
        result = drive_agent("synthetic-run", "synthetic", coverage, model, lambda: False, lambda *_: None)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(len(coverage.calls), 2)

    def test_finish_without_adjustment_and_early_stop(self):
        coverage, model = FakeCoverage(), FakeModel([decision("finish", None)])
        result = drive_agent("synthetic-run", "synthetic", coverage, model, lambda: False, lambda *_: None)
        self.assertEqual(result["status"], "completed_agent")
        self.assertEqual(len(coverage.calls), 1)
        coverage, model = FakeCoverage(), FakeModel([])
        result = drive_agent("synthetic-run", "synthetic", coverage, model, lambda: True, lambda *_: None)
        self.assertEqual(result["status"], "cancelled")
        self.assertEqual(coverage.calls, [])
        self.assertEqual(model.calls, [])


if __name__ == "__main__": unittest.main()
