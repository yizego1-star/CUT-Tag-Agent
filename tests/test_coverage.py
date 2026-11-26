"""Failure/control checks use synthetic local files, not the private BAM."""
import json
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from cuttag_web.coverage import DockerCoverageRunner
from cuttag_web.server import Workbench


class CoverageTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        sample = self.root / "data/private/synthetic"
        (sample / "processed").mkdir(parents=True)
        (sample / "qc").mkdir()
        (sample / "provenance").mkdir()
        self.bam = sample / "processed/synthetic.hg38.dedup.bam"
        self.bam.write_bytes(b"synthetic test input, not a BAM")
        Path(str(self.bam) + ".bai").write_bytes(b"synthetic index")
        (sample / "qc/fastp.json").write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 2}}}))
        self.runner = DockerCoverageRunner(self.root)
        self.workbench = Workbench(self.root, coverage_runner=self.runner)

    def tearDown(self):
        self.temp.cleanup()

    def wait_done(self, run_id):
        until = time.monotonic() + 3
        while time.monotonic() < until:
            item = self.workbench.get_run(run_id)
            if item["status"] not in ("queued", "running"):
                return item
            time.sleep(.01)
        self.fail("Synthetic control did not finish")

    def test_fixed_commands_and_readonly_input_mount(self):
        folder = self.runner.output_dir("run-123456789abc")
        command = self.runner.command("run-123456789abc", self.bam, folder, "sha256:synthetic")
        self.assertIn("type=bind,src=" + str(self.bam.parent) + ",dst=/input,readonly", command)
        self.assertEqual(command[command.index("--network") + 1], "none")
        self.assertEqual(command[command.index("--memory") + 1], "2g")
        self.assertIn("--read-only", command)
        self.assertNotIn(".env", " ".join(command))
        for invalid in ("../escape", "run-a;touch-file", "run-123456789abc/other"):
            with self.assertRaises(ValueError):
                self.runner.output_dir(invalid)

    def test_refuses_extra_real_parameters_and_outside_symlinks(self):
        for payload in ({"region": "chr2"}, {"target_bin_size": 50}, {"command": "anything"}):
            with self.assertRaises(ValueError):
                self.workbench.create_run({"mode": "coverage", "sample_id": "synthetic", **payload})
        index = Path(str(self.bam) + ".bai")
        index.unlink()
        (self.root / "outside.bai").write_bytes(b"synthetic")
        index.symlink_to(self.root / "outside.bai")
        with self.assertRaises(ValueError):
            self.runner.input_bam("synthetic")

    def test_timeout_stops_owned_container_and_does_not_publish_output(self):
        self.runner.timeout = 0
        process = Mock(returncode=None)
        process.poll.side_effect = lambda: process.returncode
        process.wait.side_effect = lambda **_: process.returncode
        def stop(_):
            process.returncode = 137
        with patch.object(self.runner, "probe", return_value={"ready": True, "image_id": "synthetic"}), \
             patch("cuttag_web.coverage.subprocess.Popen", return_value=process), \
             patch.object(self.runner, "cancel", side_effect=stop) as cancel:
            result = self.runner.run("run-123456789abc", "synthetic", lambda: False, lambda *_: None)
        self.assertEqual(result["status"], "timed_out")
        cancel.assert_any_call("run-123456789abc")
        self.assertEqual(result["actual_tool_executions"], 0)
        self.assertNotIn("result", result)
        self.assertTrue(result["execution"]["inputs_unchanged"])

    def test_failure_state_and_artifact_whitelist(self):
        def synthetic_failure(run_id, sample, cancelled, phase):
            phase("inputs", {})
            raise ValueError("Synthetic input failure")
        with patch.object(self.runner, "probe", return_value={"ready": True}), \
             patch.object(self.runner, "run", side_effect=synthetic_failure):
            run = self.workbench.create_run({"sample_id": "synthetic", "mode": "coverage"})
            result = self.wait_done(run["id"])
        self.assertEqual(result["status"], "failed")
        self.assertIn("Synthetic", result["error"])
        self.assertEqual(result["actual_tool_executions"], 0)
        for name in ("coverage.bw", "../.env", "source.bam", "tool-start.json"):
            with self.assertRaises(ValueError):
                self.workbench.artifact(run["id"], name)

    def test_cancellation_and_single_real_job_limit(self):
        entered = threading.Event()
        def synthetic_wait(run_id, sample, cancelled, phase):
            phase("coverage", {})
            entered.set()
            until = time.monotonic() + 2
            while not cancelled() and time.monotonic() < until:
                time.sleep(.01)
            return {"status": "cancelled", "actual_tool_executions": 0}
        with patch.object(self.runner, "probe", return_value={"ready": True}), \
             patch.object(self.runner, "run", side_effect=synthetic_wait):
            run = self.workbench.create_run({"sample_id": "synthetic", "mode": "coverage"})
            self.assertTrue(entered.wait(1))
            with self.assertRaises(ValueError):
                self.workbench.create_run({"sample_id": "synthetic", "mode": "coverage"})
            self.workbench.cancel_run(run["id"])
            result = self.wait_done(run["id"])
        self.assertEqual(result["status"], "cancelled")
        self.assertNotIn("result", result)

    def test_incorrect_result_cannot_be_marked_successful(self):
        folder = self.root / "bad-output"
        folder.mkdir()
        (folder / "result.json").write_text(json.dumps({"validation": "passed", "chrom": "chr2", "start": 0, "end": 100, "bin_size": 25}))
        with self.assertRaises(ValueError):
            self.runner.validate_result(folder)


if __name__ == "__main__":
    unittest.main()
