"""Run a bounded normalization or IgG preset with a previously verified local review."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import uuid

from cuttag_web.coverage import read_record
from cuttag_web.model import ModelClient
from cuttag_web.pipeline import PipelineRunner
from tools.spikein_track import matched_inputs, run_tracks


def review_source(root, review_id):
    if not isinstance(review_id, str) or not re.fullmatch(r'run-[0-9a-f]{12}', review_id):
        raise ValueError('复核运行号格式无效')
    review = read_record(root / 'outputs/web_runs' / (review_id + '.json'))
    if review.get('id') != review_id or review.get('mode') != 'review' or review.get('status') != 'completed_review':
        raise ValueError('请选择已完成的本地质控复核')
    summaries = review.get('track_review', {}).get('sent_summaries', [])
    if not summaries:
        raise ValueError('该复核缺少已发送的质控摘要')
    return review, summaries[0]


def run_preset(root, preset, review_id, *, run_id=None, model=None, cancelled=lambda: False, update=lambda _: None):
    review, payload = review_source(root, review_id)
    sample, baseline = review['sample_id'], review['baseline_id']
    folder = root / 'outputs/preset_runs' / (run_id or 'run-' + uuid.uuid4().hex[:12])
    folder.mkdir(parents=True)
    record = {'id': folder.name, 'preset': preset, 'review_id': review_id, 'baseline_id': baseline,
              'sample_id': sample, 'model_requests': 0, 'peak_executions': 0, 'fastq_executions': 0,
              'created_at': datetime.now(timezone.utc).isoformat(), 'status': 'failed'}

    def sent():
        if cancelled():
            raise ValueError('任务已停止')
        record.update(model_requests=1, sent_summary=payload)
        update(record)

    try:
        if preset == 'normalization':
            inputs = matched_inputs(root, sample)
            spikein = payload['qc'].get('spikein', {})
            if inputs[4] != spikein.get('spikein_fragments') or inputs[6] != spikein.get('human_fragments'):
                raise ValueError('本地校验计数与质控摘要不一致')
        else:
            if payload['qc']['experiment_context']['matched_igg'] != 'unavailable':
                raise ValueError('此预设用于核对已确认无匹配 IgG 的样本')
            baseline_folder = PipelineRunner(root).output_dir(baseline)
            source = read_record(baseline_folder / 'summary.json').get('source')
            if source not in ('resume', 'fastq', 'bam'):
                raise ValueError('找不到可核对的原始分析结果')
            _, artifacts = PipelineRunner.validate(baseline_folder, source)
            record['preserved_peaks'] = [item for item in artifacts if item['name'] in ('narrow.peaks.bed', 'broad.peaks.bed')]
            if len(record['preserved_peaks']) != 2:
                raise ValueError('已有候选峰文件缺失')

        answer = (model or ModelClient(root)).review_preset(payload, sample, preset, on_request=sent)
        record['answer'] = answer
        action = answer['decision']['action']
        update(record)
        if cancelled():
            raise ValueError('任务已停止')
        if preset == 'normalization' and action == 'generate_pair':
            tracks = run_tracks(sample, folder / 'tracks', inputs, compare=True, cancelled=cancelled)
            if cancelled():
                raise ValueError('任务已停止')
            if tracks['status'] != 'completed_comparison':
                raise ValueError(tracks['error'])
            record.update(status='completed_comparison', coverage_executions=2,
                          track_record='tracks/record.json', execution_seconds=tracks['elapsed_seconds'],
                          comparison_ratio=tracks['result']['local_spikein_to_cpm_ratio'],
                          compared_preview_bins=tracks['result']['compared_preview_bins'])
        elif preset == 'igg' and action == 'keep_existing_peaks':
            record.update(status='blocked_no_matched_igg', coverage_executions=0)
        else:
            raise ValueError('模型没有给出此预设的可执行建议')
    except (OSError, ValueError, KeyError, IndexError, TypeError, KeyboardInterrupt) as exc:
        record['error'] = str(exc) or '用户停止'
        if cancelled() or isinstance(exc, KeyboardInterrupt):
            record['status'] = 'cancelled'
    finally:
        record['updated_at'] = datetime.now(timezone.utc).isoformat()
        path = folder / 'record.json'
        temporary = path.with_suffix('.tmp')
        temporary.write_text(json.dumps(record, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
        temporary.replace(path)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('preset', choices=('normalization', 'igg'))
    parser.add_argument('--review', required=True, help='已完成的本地质控复核运行号')
    args = parser.parse_args()
    if not re.fullmatch(r'run-[0-9a-f]{12}', args.review):
        parser.error('运行号格式无效')
    try:
        root = Path(__file__).resolve().parents[1]
        run_id = 'run-' + uuid.uuid4().hex[:12]
        print(root / 'outputs/preset_runs' / run_id, flush=True)
        record = run_preset(root, args.preset, args.review, run_id=run_id)
        print(record.get('error', record['status']), flush=True)
        return int(record['status'] in ('failed', 'cancelled'))
    except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
