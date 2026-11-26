"""Terminal handoff for one existing local run. Never starts analysis or calls a model."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
BACKGROUND_LABEL = None


def read_run(run_id):
    if not re.fullmatch(r'run-[0-9a-f]{12}', run_id):
        raise ValueError('运行编号无效')
    path = (ROOT / 'outputs/web_runs' / (run_id + '.json')).resolve()
    if path.parent != (ROOT / 'outputs/web_runs').resolve() or path.stat().st_size > 2 * 1024 ** 2:
        raise ValueError('运行记录路径或大小无效')
    return json.loads(path.read_text())


def emit(message):
    print(datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S') + '  ' + message, flush=True)


def main():
    global BACKGROUND_LABEL
    parser = argparse.ArgumentParser()
    parser.add_argument('run_id')
    parser.add_argument('--status', action='store_true')
    args = parser.parse_args()
    if not args.status and re.fullmatch(r'run-[0-9a-f]{12}', args.run_id):
        BACKGROUND_LABEL = 'local.cuttag.watch.' + args.run_id
    last, inactive_checks = None, 0
    while True:
        item = read_run(args.run_id)
        state = (item['status'], item.get('message', ''))
        if state != last:
            emit(item['status'] + ' | ' + item.get('message', ''))
            last = state
        if item['status'] not in ('queued', 'running'):
            if item['status'] == 'completed_pipeline':
                emit('基础分析完成。结果在 outputs/pipeline_runs/' + args.run_id)
                emit('此状态表示工具与文件检查完成，不代表实验质量合格。')
            else:
                emit('任务未成功完成，请保留日志并告知编码助手。' + item.get('error', ''))
            return 0 if item['status'] == 'completed_pipeline' else 1
        if args.status:
            return 0
        # A server that vanished must not leave an eternal stale success-looking watcher.
        try:
            result = subprocess.run(['lsof', '-nP', '-iTCP:8766', '-sTCP:LISTEN', '-t'],
                                    capture_output=True, text=True, timeout=5)
            inactive_checks = 0 if result.returncode == 0 and result.stdout.strip() else inactive_checks + 1
        except (OSError, subprocess.TimeoutExpired):
            inactive_checks += 1
        if inactive_checks >= 3:
            emit('本地服务没有运行；本脚本不会重复启动计算。请保留当前 Docker 状态并告知编码助手。')
            return 2
        time.sleep(5)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError) as exc:
        emit('无法读取本地运行状态：' + str(exc))
        sys.exit(2)

    finally:
        # launchctl submit may retain a submitted job; remove it on every terminal outcome.
        # Removing our own registration also releases the parent caffeinate assertions.
        if BACKGROUND_LABEL:
            try:
                subprocess.run(['/bin/launchctl', 'remove', BACKGROUND_LABEL],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=5)
            except (OSError, subprocess.TimeoutExpired):
                pass
