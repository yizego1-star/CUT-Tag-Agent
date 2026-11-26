#!/bin/bash
# Attach a macOS-owned keep-awake/status task to an existing run. Never starts analysis.
set -euo pipefail
cuttag_root="$(cd "$(dirname "$0")/.." && pwd)"
cuttag_run="${1:?用法：bash tools/terminal_handoff.sh <run-id>}"
if [[ ! "$cuttag_run" =~ ^run-[0-9a-f]{12}$ ]]; then
  echo '运行编号无效' >&2
  exit 2
fi
cd "$cuttag_root"
/usr/bin/python3 tools/watch_local_run.py "$cuttag_run" --status
mkdir -p .local/terminal_watch outputs/terminal_runs
cuttag_label="local.cuttag.watch.$cuttag_run"
cuttag_log="$cuttag_root/outputs/terminal_runs/$cuttag_run.log"
if /bin/launchctl list "$cuttag_label" >/dev/null 2>&1; then
  echo '本次后台保护已经登记；不会重复启动分析。'
  echo "查看进度：tail -f '$cuttag_log'"
  exit 0
fi
/bin/launchctl submit -l "$cuttag_label" -o "$cuttag_log" -e "$cuttag_log" -- /usr/bin/caffeinate -ims /usr/bin/python3 "$cuttag_root/tools/watch_local_run.py" "$cuttag_run"
printf '%s\n' "$cuttag_label" > "$cuttag_root/.local/terminal_watch/$cuttag_run.label"
echo '已交 macOS 后台任务服务：防闲置休眠并记录进度，不重复运行分析。'
echo '完成后自动解除本脚本的防休眠。请保持 Docker 和当前本地服务运行。'
echo '合盖工作需保持电源、外屏和外接键鼠；防闲置休眠不能代替这些条件。'
echo "查看进度：tail -f '$cuttag_log'"
