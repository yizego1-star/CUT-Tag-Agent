"""One model decision on a completed local CUT&Tag coverage track. No full rerun."""
import argparse
import json
from pathlib import Path
import uuid

from cuttag_web.coverage import read_record
from cuttag_web.model import ModelClient
from cuttag_web.pipeline import PipelineRunner
from cuttag_web.track_review import review_baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline_run')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    folder = PipelineRunner(root).output_dir(args.baseline_run)
    finished = read_record(folder / 'resume-result.json')
    if finished.get('status') != 'completed_pipeline' or finished.get('execution', {}).get('exit_code') != 0:
        parser.exit(1, '没有通过检查的本地基础运行，未调用模型。\n')
    parent = finished.get('execution', {}).get('resumed_from', '')
    PipelineRunner.container_name(parent)
    original = read_record(root / 'outputs/web_runs' / (parent + '.json'))
    scope = read_record(root / '.local/model_summary_authorization.json')
    if original.get('sample_id') not in scope.get('sample_ids', []):
        parser.exit(1, '样本未列入最小摘要的许可范围，未调用模型。\n')
    result_dir = root / 'outputs' / 'agent_reviews'
    result_dir.mkdir(exist_ok=True)
    receipt = result_dir / ('review-' + uuid.uuid4().hex[:12] + '.json')
    result = review_baseline(folder, 'resume', ModelClient(root, timeout=30))
    result['baseline_run'] = args.baseline_run
    receipt.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2) + '\n')
    print(receipt)
    return int(result['status'] != 'completed_review')


if __name__ == '__main__':
    raise SystemExit(main())
