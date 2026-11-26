"""Prepare a project-local H9 tool environment for Linux."""
import argparse
import ctypes
import hashlib
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile


ROOT = Path(__file__).resolve().parent.parent
PIXI_URL = 'https://github.com/prefix-dev/pixi/releases/download/v0.77.1/pixi-x86_64-unknown-linux-musl.tar.gz'
PIXI_SHA256 = '74dbe15255c763396ac847eb455573c6400e9caf06614c3fc00cb3de6d1099d5'
PIXI_BINARY_SHA256 = '5115a89a9189a2e4e7e8d2f04236a7be586d8f6091dfc9ea869fb3c4a52b6935'
MANIFEST = """[workspace]
name = "cuttag-h9-linux"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[dependencies]
python = "3.11.*"
fastp = "0.24.0.*"
bowtie2 = "2.5.4.*"
samtools = "1.21.*"
macs3 = "3.0.2.*"
deeptools = "3.5.6.*"
pysam = "*"
pybigwig = "*"
"""


class PathRule(ctypes.Structure):
    _pack_ = 1
    _fields_ = [('access', ctypes.c_uint64), ('parent_fd', ctypes.c_int32)]


def restrict_writes(root):
    """Permit filesystem writes, creates, deletes and renames only below root."""
    libc = ctypes.CDLL(None, use_errno=True)

    def check(result):
        if result < 0:
            code = ctypes.get_errno()
            raise OSError(code, os.strerror(code))
        return result

    # x86_64 Linux uses the generic Landlock syscall numbers.
    create, add, restrict = 444, 445, 446
    try:
        abi = check(libc.syscall(create, 0, 0, 1))
        if abi < 3:
            raise ValueError(f'Landlock ABI {abi} 不能限制文件截断')
        # WRITE_FILE, REMOVE_*, MAKE_*, REFER and TRUNCATE.
        mask = (1 << 1) | sum(1 << bit for bit in range(4, 15))
        rights = ctypes.c_uint64(mask)
        ruleset = check(libc.syscall(create, ctypes.byref(rights), ctypes.sizeof(rights), 0))
        try:
            fd = os.open(root, os.O_PATH | os.O_CLOEXEC | os.O_DIRECTORY)
            try:
                rule = PathRule(mask, fd)
                check(libc.syscall(add, ruleset, 1, ctypes.byref(rule), 0))
                check(libc.prctl(38, 1, 0, 0, 0))  # PR_SET_NO_NEW_PRIVS
                check(libc.syscall(restrict, ruleset, 0))
            finally:
                os.close(fd)
        finally:
            os.close(ruleset)
    except OSError as exc:
        raise ValueError(f'目录写入限制未启用，已停止。系统返回 {exc.strerror}') from exc
    return abi


def ensure_text(path, expected):
    if path.is_symlink() or (path.exists() and path.read_text() != expected):
        raise ValueError(f'已有不同的配置，保留原文件 {path}')
    if not path.exists():
        with path.open('x') as stream:
            stream.write(expected)


def install_pixi(binary, scratch, env):
    if binary.exists():
        if binary.is_symlink() or hashlib.sha256(binary.read_bytes()).hexdigest() != PIXI_BINARY_SHA256:
            raise ValueError(f'已有 Pixi 文件校验失败，保留 {binary}')
        return
    with tempfile.TemporaryDirectory(dir=scratch) as staging:
        staging = Path(staging)
        archive = staging / 'pixi.tar.gz'
        subprocess.run(['/usr/bin/curl', '-q', '-fLSs', '--proto', '=https', '--proto-redir', '=https',
                        '--connect-timeout', '15', '--max-time', '120', '--retry', '2', '-o', str(archive), PIXI_URL],
                       env=env, check=True)
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        if digest != PIXI_SHA256:
            raise ValueError('Pixi 下载 SHA-256 不匹配')
        candidate = staging / 'pixi'
        with tarfile.open(archive, 'r:gz') as bundle:
            members = [m for m in bundle.getmembers() if m.isfile() and Path(m.name).name == 'pixi']
            if len(members) != 1:
                raise ValueError('Pixi 下载包的文件清单无效')
            with bundle.extractfile(members[0]) as source, candidate.open('xb') as target:
                shutil.copyfileobj(source, target)
        if hashlib.sha256(candidate.read_bytes()).hexdigest() != PIXI_BINARY_SHA256:
            raise ValueError('Pixi 可执行文件 SHA-256 不匹配')
        candidate.chmod(0o700)
        candidate.replace(binary)


def prepare_runtime(root):
    runtime = root / '.h9-linux'
    scratch = runtime / 'tmp'
    workspace = runtime / 'workspace'
    config = workspace / '.pixi'
    for directory in (scratch, config, runtime / 'bin'):
        if directory.resolve() != directory:
            raise ValueError(f'安装路径不能经过符号链接 {directory}')
        directory.mkdir(mode=0o700, parents=True, exist_ok=True)
    network_settings = {'http_proxy', 'https_proxy', 'all_proxy', 'no_proxy', 'ssl_cert_file', 'ssl_cert_dir'}
    env = {key: value for key, value in os.environ.items() if key.lower() in network_settings}
    env.update(PATH='/usr/bin:/bin', TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch),
               PIXI_HOME=str(runtime / 'pixi-home'), PIXI_CACHE_DIR=str(runtime / 'cache'),
               PIXI_NO_CONFIG='1', PIXI_NO_TELEMETRY='1', RATTLER_AUTH_FILE=str(runtime / 'auth.json'),
               XDG_CACHE_HOME=str(runtime / 'xdg-cache'), XDG_CONFIG_HOME=str(runtime / 'xdg-config'),
               XDG_DATA_HOME=str(runtime / 'xdg-data'), XDG_STATE_HOME=str(runtime / 'xdg-state'),
               MPLCONFIGDIR=str(runtime / 'matplotlib'), PYTHONDONTWRITEBYTECODE='1')
    ensure_text(runtime / 'auth.json', '{}\n')
    ensure_text(config / 'config.toml', '[cache]\nnetfs-redirect = "never"\n')
    manifest = workspace / 'pixi.toml'
    ensure_text(manifest, MANIFEST)
    pixi = runtime / 'bin/pixi'
    install_pixi(pixi, scratch, env)
    return workspace, pixi, env


def install(root):
    workspace, pixi, env = prepare_runtime(root)
    if subprocess.check_output([str(pixi), '--version'], env=env, text=True).split()[:2] != ['pixi', '0.77.1']:
        raise ValueError(f'已存在的 Pixi 版本不符，保留 {pixi}')
    subprocess.run([str(pixi), 'install', '--manifest-path', str(workspace / 'pixi.toml')], env=env,
                   cwd=workspace, check=True)
    binaries = workspace / '.pixi/envs/default/bin'
    env['PATH'] = str(binaries) + ':/usr/bin:/bin'
    for name in ('fastp', 'bowtie2', 'samtools', 'macs3', 'bamCoverage'):
        result = subprocess.check_output([str(binaries / name), '--version'], env=env, text=True,
                                         stderr=subprocess.STDOUT)
        print(result.splitlines()[0], flush=True)
    subprocess.run([str(binaries / 'python'), '-B', '-c', 'import pysam, pyBigWig; print("Python imports OK")'],
                   env=env, check=True)
    print('H9 工具环境可用', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--probe', action='store_true', help='只检查项目目录能否创建和移动文件')
    modes.add_argument('--prepare', action='store_true', help='下载并校验 Pixi，写入项目配置，不运行包安装器')
    args = parser.parse_args()
    if platform.system() != 'Linux' or platform.machine() != 'x86_64':
        parser.error('此脚本仅适用于 Linux x86_64')
    if Path.cwd() != ROOT or ROOT.resolve() != ROOT:
        parser.error(f'请进入项目目录 {ROOT}。路径不能经过符号链接')
    os.umask(0o077)
    try:
        # The full-install mode launches a package manager. The probe and preparation
        # modes only write their explicit project-local paths.
        if not (args.probe or args.prepare):
            abi = restrict_writes(ROOT)
            print(f'包安装进程已启用 Landlock ABI {abi}', flush=True)
        with tempfile.TemporaryDirectory(prefix='.h9-write-check-', dir=ROOT) as check:
            sample = Path(check) / 'before'
            sample.write_text('ok')
            sample.rename(sample.with_name('after'))
        print('项目目录写入检查通过', flush=True)
        if args.prepare:
            _, pixi, _ = prepare_runtime(ROOT)
            print(f'Pixi 已下载并校验：{pixi}', flush=True)
        elif not args.probe:
            install(ROOT)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        sys.exit(str(exc))


if __name__ == '__main__':
    main()
