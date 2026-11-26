"""Run the fixed downstream pipeline from a recorded local FASTQ checkpoint."""
import argparse
import fcntl
import json
from pathlib import Path
import signal
import sys
import uuid

from cuttag_web.pipeline import PipelineRunner, RESUME_IMAGE


def main():
    parser = argparse.ArgumentParser(description='从已停止的本地 FASTQ 任务的 BAM 续跑。')
    parser.add_argument('previous_run')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runner = PipelineRunner(root, image=RESUME_IMAGE)
    sample_id, _ = runner.resume_inputs(args.previous_run)
    run_id = 'run-' + uuid.uuid4().hex[:12]
    folder = runner.output_dir(run_id)
    stop = False

    def request_stop(*_):
        nonlocal stop
        stop = True

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)

    def phase(_key, message, _execution):
        print(message, flush=True)

    with (folder.parent / '.resume.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('已有本地续跑任务正在运行') from None
        print(f'续跑来源：{args.previous_run}\n本轮结果目录：{folder}\n工具日志：{folder / "tool.log"}', flush=True)
        result = runner.run(run_id, sample_id, 'resume', lambda: stop, phase, resume_from=args.previous_run)
        receipt = folder / 'resume-result.json'
        receipt.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
        if result['status'] == 'completed_pipeline':
            print(f'验收通过。结果说明：{folder / "report.md"}', flush=True)
        else:
            print(f"未完成：{result['execution'].get('error', result['status'])}。请查看 {receipt}", file=sys.stderr)
        if result['execution'].get('scratch_cleanup_error'):
            print(f'临时目录清理失败，请检查 {folder / ".scratch"}', file=sys.stderr)
        return 0 if result['status'] == 'completed_pipeline' else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        sys.exit(str(exc))
