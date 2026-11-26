#!/usr/bin/env bash
# Launch the public H9 workflow from a Linux checkout.
set -euo pipefail
export PATH=/usr/bin:/bin
umask 077
root="$(cd -- "$(dirname -- "$0")/.." && pwd -P)"
runtime="$root/.h9-linux"
state="$runtime/h9-session"
logs="$root/outputs/h9_development/native-logs"
python="$runtime/offline-tools/.pixi/envs/default/bin/python"
socket="$state/h9.sock"

cd -- "$root"
[[ "$(pwd -P)" == "$root" && -x "$python" ]] || { echo '请先在项目目录准备离线工具环境。' >&2; exit 1; }
for path in "$root/tools" "$state" "$state/home" "$state/tmp" "$logs"; do
    [[ "$(realpath -m "$path")" == "$path" ]] || { echo "路径经过符号链接：$path" >&2; exit 1; }
done
mkdir -p -- "$state/home" "$state/tmp" "$logs"
case "${1:-}" in
  smoke) flag=--prepare-and-smoke ;;
  full) flag=--full ;;
  *) echo '用法：bash tools/h9-job.sh smoke 或 full' >&2; exit 2 ;;
esac
mode="$1"
log="$logs/$mode.log"
[[ ! -L "$log" && ! -L "$socket" ]] || { echo '日志或套接字为符号链接。' >&2; exit 1; }
if [[ "$mode" == full && ! -f "$root/outputs/h9_development/prefix-100000/COMPLETE.json" ]]; then
    echo '请先等待限量试跑完成。' >&2
    exit 1
fi

if [[ "${2:-}" == --foreground ]]; then
    cd -- "$root"
    exec >> "$log" 2>&1
    printf '\n[%s] H9 %s 开始\n' "$(date -Is)" "$mode"
    status=0
    env -i PATH=/usr/bin:/bin HOME="$state/home" TMPDIR="$state/tmp" \
        nice -n 10 "$python" -E -s -B -m tools.h9_native "$flag" || status=$?
    printf '[%s] H9 %s 退出码 %s\n' "$(date -Is)" "$mode" "$status"
    exit "$status"
fi

if command -v tmux >/dev/null 2>&1; then
    env -i PATH=/usr/bin:/bin HOME="$state/home" TMPDIR="$state/tmp" SHELL=/bin/bash TERM=xterm \
        tmux -f /dev/null -S "$socket" new-session -d -s "h9-$mode" -c "$root" \
        /bin/bash "$root/tools/h9-job.sh" "$mode" --foreground
    printf '已启动 tmux 会话 h9-%s\n' "$mode"
    printf '进入会话：tmux -S %q attach -t h9-%s\n' "$socket" "$mode"
else
    nohup /bin/bash "$root/tools/h9-job.sh" "$mode" --foreground </dev/null >> "$log" 2>&1 &
    printf '已通过 nohup 启动，PID %s\n' "$!"
fi
printf '日志：%s\n' "$log"
