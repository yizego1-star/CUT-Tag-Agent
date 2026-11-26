"""Local API checks with synthetic reports; never change the development sample."""
import http.client
import json
import tempfile
import threading
import time
import unittest
from http.server import ThreadingHTTPServer
from pathlib import Path

from cuttag_web.server import PROJECT_ROOT, Workbench, handler_for


class WorkbenchTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.sample = self.root / "data" / "private" / "example"
        (self.sample / "qc").mkdir(parents=True)
        (self.sample / "provenance").mkdir()
        self.write_json("qc/fastp.json", {"summary": {
            "before_filtering": {"total_reads": 100, "q30_rate": .91},
            "after_filtering": {"total_reads": 96},
        }})
        self.write_json("provenance/source_manifest.json", {"files": [{
            "local_relative_path": "qc/fastp.json",
            "expected_bytes": (self.sample / "qc" / "fastp.json").stat().st_size,
        }]})
        self.write_json("provenance/verification.json", {"status": "passed", "checks": [{"name": "example"}]})
        (self.sample / "qc" / "hg38.dedup.fragment_length_counts.txt").write_text("3 25\n7 152\n")
        (self.root / ".env").write_text("CUTTAG_LLM_API_KEY=synthetic-secret-for-test\nCUTTAG_LLM_MODEL=example-model\n")
        self.before = self.input_snapshot()
        self.workbench = Workbench(self.root, step_delay=.025)
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), handler_for(self.workbench, PROJECT_ROOT / "web" / "dist"))
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def tearDown(self):
        self.assertEqual(self.before, self.input_snapshot(), "Input files must remain unchanged")
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.directory.cleanup()

    def write_json(self, name, data):
        (self.sample / name).write_text(json.dumps(data))

    def input_snapshot(self):
        return {str(path.relative_to(self.sample)): (path.read_bytes(), path.stat().st_mtime_ns)
                for path in self.sample.rglob("*") if path.is_file()}

    def request(self, path, payload=None, origin=None):
        connection = http.client.HTTPConnection("127.0.0.1", self.server.server_port, timeout=3)
        headers = {"Content-Type": "application/json"}
        if origin:
            headers["Origin"] = origin
        connection.request("POST" if payload is not None else "GET", path,
                           json.dumps(payload) if payload is not None else None, headers)
        response = connection.getresponse()
        content = response.read()
        result = response.status, dict(response.getheaders()), content
        connection.close()
        return result

    def wait_finished(self, run_id):
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            item = self.workbench.get_run(run_id)
            if item["status"] not in ("queued", "running"):
                return item
            time.sleep(.01)
        self.fail("Demo did not finish")

    def test_reports_and_configuration_do_not_expose_secrets(self):
        status, _, body = self.request("/api/samples")
        self.assertEqual(status, 200)
        sample = json.loads(body)["samples"][0]
        self.assertEqual(sample["metrics"]["total_reads"], 100)
        self.assertEqual(sample["fragment_distribution"]["counts"][2], 3)
        self.assertEqual(sample["fragment_distribution"]["counts"][15], 7)
        _, _, body = self.request("/api/health")
        self.assertTrue(json.loads(body)["model_key_configured"])
        self.assertNotIn(b"synthetic-secret", body)
        self.assertFalse(json.loads(body)["model_requests_sent"])

    def test_private_files_are_not_web_assets(self):
        for path in ("/.env", "/data/private/example/qc/fastp.json", "/../AGENTS.md", "/%2e%2e/.env"):
            with self.subTest(path=path):
                self.assertEqual(self.request(path)[0], 404)

    def test_rejects_real_execution_and_cross_origin_creation(self):
        self.assertEqual(self.request("/api/runs", {"sample_id": "example", "mode": "real"})[0], 400)
        self.assertEqual(self.request("/api/runs", {"sample_id": "example", "command": "anything"})[0], 400)
        self.assertEqual(self.request("/api/runs", {"sample_id": "example", "target_bin_size": True})[0], 400)
        self.assertEqual(self.request("/api/runs", {"sample_id": "example"}, "https://example.com")[0], 403)
        self.assertEqual(self.workbench.list_runs(), [])

    def test_demo_records_plan_without_fabricating_analysis(self):
        status, _, body = self.request("/api/runs", {"sample_id": "example", "name": "Example demo", "target_bin_size": 50})
        self.assertEqual(status, 201)
        run_id = json.loads(body)["id"]
        item = self.wait_finished(run_id)
        self.assertEqual(item["status"], "completed_demo")
        self.assertEqual(item["actual_tool_executions"], 0)
        self.assertEqual(item["actual_model_requests"], 0)
        self.assertEqual(item["plan"]["source"], "fixed_demo")
        self.assertEqual(item["plan"]["bin_size_after"], 50)
        self.assertNotIn("new_metrics", item)
        self.assertNotIn("output_bw", item)
        status, headers, body = self.request("/api/runs/" + run_id + "/report")
        self.assertEqual(status, 200)
        self.assertIn("attachment", headers["Content-Disposition"])
        self.assertEqual(json.loads(body), item)
        restored = Workbench(self.root)
        self.assertEqual(restored.get_run(run_id), item)

    def test_stopping_last_step_is_respected(self):
        self.workbench.step_delay = .07
        item = self.workbench.create_run({"sample_id": "example"})
        run_id = item["id"]
        deadline = time.monotonic() + 3
        while time.monotonic() < deadline:
            item = self.workbench.get_run(run_id)
            if item["steps"][-1]["status"] == "running":
                break
            time.sleep(.002)
        self.workbench.cancel_run(run_id)
        self.assertEqual(self.wait_finished(run_id)["status"], "cancelled")


if __name__ == "__main__":
    unittest.main()
