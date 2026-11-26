"""Run the public H9 workflow inside this checkout on Linux."""
import argparse
import fcntl
import gzip
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from xml.sax.saxutils import escape

from tools import run_h9
from tools.h9_samples import SAMPLES


ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT
RUNTIME = ROOT / '.h9-linux'
PREFIX = RUNTIME / 'offline-tools/.pixi/envs/default'
DATA = ROOT / 'data/public/GSE200768_H9_H3K4me1'
REFERENCE = ROOT / 'data/support/shared-reference/reference/hg38'
OUTPUT = ROOT / 'outputs/h9_development'
FASTA_MD5 = '6d3c82e1e12b127d526395294526b9c8'
INDEX_SUFFIXES = ('1', '2', '3', '4', 'rev.1', 'rev.2')


def inside(path):
    if path.resolve() != path or ROOT not in path.parents:
        raise ValueError(f'路径不在项目目录内或经过符号链接：{path}')
    return path


def isolated_env(scratch):
    inside(scratch).mkdir(parents=True, exist_ok=True)
    home, cache, config = (inside(scratch / name) for name in ('home', 'cache', 'config'))
    for path in (home, cache, config):
        path.mkdir(exist_ok=True)
    font_config = inside(scratch / 'fonts.conf')
    font_text = f'<fontconfig><cachedir>{escape(str(cache))}</cachedir></fontconfig>\n'
    if font_config.exists():
        if font_config.read_text() != font_text:
            raise ValueError('字体缓存配置不匹配')
    else:
        font_config.write_text(font_text)
    return {'PATH': str(PREFIX / 'bin') + ':/usr/bin:/bin', 'HOME': str(home), 'PWD': str(CODE),
            'TMPDIR': str(scratch), 'TMP': str(scratch), 'TEMP': str(scratch),
            'XDG_CACHE_HOME': str(cache), 'XDG_CONFIG_HOME': str(config),
            'XDG_DATA_HOME': str(scratch / 'data'), 'XDG_STATE_HOME': str(scratch / 'state'),
            'XDG_RUNTIME_DIR': str(scratch), 'MPLCONFIGDIR': str(cache), 'MPLBACKEND': 'Agg',
            'FONTCONFIG_FILE': str(font_config), 'PYTHONDONTWRITEBYTECODE': '1',
            'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'LC_ALL': 'C', 'PYTHONPATH': str(CODE)}


def hashes(path):
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with inside(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def verify_fastqs():
    manifest = run_h9.read(inside(DATA / 'metadata/raw_download_manifest.json'))
    checked = run_h9.read(inside(DATA / 'metadata/fastq_validation.source.json'))
    rows = {row['local_relative_path']: row for row in manifest['files']}
    expected = {f'raw/{sample}/{sample}_{mate}.fastq.gz' for sample in SAMPLES for mate in (1, 2)}
    if (manifest['source_is_public'] is not True or manifest['status'] != 'downloaded' or set(rows) != expected
            or checked['status'] != 'passed' or {row['file'] for row in checked['files']} != expected):
        raise ValueError('H9 公开数据清单不完整')
    state = {}
    for name in sorted(rows):
        row = rows[name]
        file = inside(DATA / name)
        before = file.stat()
        if before.st_size != row['bytes'] or hashes(file) != (row['md5'], row['sha256']):
            raise ValueError(f'上传的 FASTQ 校验失败：{name}')
        after = file.stat()
        if (after.st_size, after.st_mtime_ns) != (before.st_size, before.st_mtime_ns):
            raise ValueError(f'FASTQ 校验过程中发生变化：{name}')
        state[name] = [after.st_size, after.st_mtime_ns]
        print(f'已核对 FASTQ {name}', flush=True)
    checked['native_transfer_verification'] = {'method': 'byte-identical to the gzip/record/pair-checked public FASTQ',
                                              'size_and_mtime_ns': state}
    run_h9.write(inside(DATA / 'metadata/fastq_validation.json'), checked)


def index_inputs(index):
    files = [inside(index / f'hg38_noalt.{part}.{extension}')
             for extension in ('bt2', 'bt2l') for part in INDEX_SUFFIXES
             if (index / f'hg38_noalt.{part}.{extension}').exists()]
    if len(files) != len(INDEX_SUFFIXES) or len({path.suffix for path in files}) != 1:
        raise ValueError('Bowtie 2 索引需要一整套 .bt2 或 .bt2l 文件')
    state = {}
    for path in files:
        info = path.stat()
        if not info.st_size:
            raise ValueError(f'Bowtie 2 索引文件为空：{path}')
        state[path.name] = [info.st_size, info.st_mtime_ns]
    return state


def verify_index(index, env):
    output = subprocess.check_output([str(PREFIX / 'bin/bowtie2-inspect'), '-s', str(index / 'hg38_noalt')],
                                     cwd=str(REFERENCE), env=env, text=True, timeout=120)
    indexed = [(f[1], int(f[2])) for row in output.splitlines() if (f := row.split()) and f[0].startswith('Sequence-')]
    sizes = [(name, int(size)) for name, size in
             (line.split() for line in inside(REFERENCE / 'hg38.bam-compatible.chrom.sizes').read_text().splitlines())]
    if len(sizes) != 195 or indexed != sizes:
        raise ValueError('Bowtie 2 索引的参考序列名称、长度或顺序不吻合')
    return index_inputs(index)


def prepare_index(env):
    index = inside(REFERENCE / 'bowtie2')
    if not index.exists():
        archive = inside(REFERENCE / 'hg38.analysisSet.fa.gz')
        if hashes(archive)[0] != FASTA_MD5:
            raise ValueError('hg38 FASTA 与 UCSC 官方 MD5 不一致')
        if shutil.disk_usage(REFERENCE).free < 30 * 1024 ** 3:
            raise ValueError('准备索引前需要至少 30 GiB 可用空间')
        log = inside(REFERENCE / 'bowtie2-build.log')
        print(f'开始构建 hg38 索引，日志：{log}', flush=True)
        with tempfile.TemporaryDirectory(prefix='.h9-index-', dir=REFERENCE) as tmp:
            staged = Path(tmp) / 'bowtie2'
            staged.mkdir()
            fasta = Path(tmp) / 'genome.fa'
            with gzip.open(archive, 'rb') as inp, fasta.open('xb') as out:
                shutil.copyfileobj(inp, out)
            with log.open('ab') as stream:
                subprocess.run([str(PREFIX / 'bin/bowtie2-build'), '--threads', '4', str(fasta), str(staged / 'hg38_noalt')],
                               cwd=tmp, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
            verify_index(staged, env)
            staged.rename(index)
    state = verify_index(index, env)
    run_h9.write(inside(REFERENCE / 'native_index_validation.json'), {'inputs': state, 'reference_sequences': 195})
    print('hg38 索引检查通过（195 条参考序列）', flush=True)


def run_scope(full):
    code_sources = ('tools/h9_native.py', 'tools/run_h9.py', 'tools/h9_worker.py',
                    'tools/pipeline/worker.py', 'tools/h9_samples.py', 'cuttag_web/h9_review.py')
    runtime_digest = inside(PREFIX.parents[2] / '.offline-bundle.sha256').read_text().strip()
    index_status = run_h9.read(inside(REFERENCE / 'native_index_validation.json'))
    if index_status['inputs'] != index_inputs(REFERENCE / 'bowtie2') or index_status['reference_sequences'] != 195:
        raise ValueError('索引已变化，请重新运行 --prepare-inputs')
    checks = run_h9.read(inside(DATA / 'metadata/fastq_validation.json'))['native_transfer_verification']
    for name, saved in checks['size_and_mtime_ns'].items():
        stat = inside(DATA / name).stat()
        if [stat.st_size, stat.st_mtime_ns] != saved:
            raise ValueError('H9 FASTQ 已变化，请重新运行 --prepare-inputs')
    recipe = {'runtime_sha256': runtime_digest, 'code': {name: hashes(CODE / name)[1] for name in code_sources},
              'index': index_status['inputs'], 'chromosome_sizes': hashes(REFERENCE / 'hg38.bam-compatible.chrom.sizes')[1]}
    gate = inside(OUTPUT / 'prefix-100000/COMPLETE.json')
    if full and (not gate.is_file() or run_h9.read(gate).get('recipe') != recipe):
        raise ValueError('请先用当前工具、代码和索引完成 --smoke')

    run_h9.ROOT, run_h9.DATA, run_h9.OUTPUT = ROOT, DATA, OUTPUT
    run_h9.REFERENCE, run_h9.INDEX = REFERENCE, REFERENCE / 'bowtie2'
    run_h9.SIZE_TABLE = REFERENCE / 'hg38.bam-compatible.chrom.sizes'
    scope = inside(OUTPUT / ('full' if full else 'prefix-100000'))
    for path in scope.rglob('*'):
        if path.is_symlink():
            raise ValueError(f'输出目录含符号链接：{path}')
    completed = run_h9.run(None if full else 100000, compare_tools=True, native_recipe=recipe)
    run_h9.write(inside(completed / 'COMPLETE.json'), {'status': 'completed', 'recipe': recipe,
                                                   'scope': 'full' if full else 'prefix_smoke'})
    print(f'计算完成：{completed / "COMPLETE.json"}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--prepare-inputs', action='store_true')
    modes.add_argument('--prepare-and-smoke', action='store_true')
    modes.add_argument('--smoke', action='store_true')
    modes.add_argument('--full', action='store_true')
    args = parser.parse_args()
    if ROOT.resolve() != ROOT or Path(__file__).resolve().parent.parent != CODE or Path.cwd() != CODE:
        raise ValueError(f'请在项目目录 {CODE} 执行')
    if Path(sys.executable).resolve() != (PREFIX / 'bin/python').resolve():
        raise ValueError('请用项目内离线环境里的 Python 执行')
    os.umask(0o077)
    inside(RUNTIME).mkdir(exist_ok=True)
    with inside(RUNTIME / 'h9-native.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.prepare_inputs or args.prepare_and_smoke:
            verify_fastqs()
            prepare_index(isolated_env(RUNTIME / 'h9-prepare-tmp'))
        if not args.prepare_inputs:
            run_scope(args.full)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        sys.exit(str(exc))
