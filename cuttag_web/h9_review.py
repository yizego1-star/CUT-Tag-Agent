"""Model contract for the verified public H9 development group."""
import math


CONTEXT = {
    'series': 'GSE200768',
    'assay': 'CUT&Tag',
    'cell_line': 'H9',
    'cell_type': 'human embryonic stem cells',
    'reference': 'hg38',
    'target': 'H3K4me1',
    'target_antibody': 'Abcam ab8895',
    'igg_antibody': 'Santa Cruz sc-2027',
    'target_replicates': 2,
    'igg_replicates': 2,
    'igg_design': 'same study and cell line; used as a shared control group, no individual pairing assumed',
    'peak_mode': 'narrow',
    'coverage_normalization': 'CPM',
    'quantitative_spikein_design': 'unconfirmed',
}
SYSTEM = """你在复核公开 H9 CUT&Tag 的候选窄峰。任务是决定是否有依据继续比较找峰，或保留结果供 IGV 人工复核。
所给 H3K4me1、抗体、H9 背景、两个目标重复和两个同研究 IgG 对照都是已知信息。
CPM 是浏览用轨道的当前标准化；定量 spike-in 设计未知，本次不做跨条件定量。不要因 spike-in 未知而单独阻止窄峰比较。
prefix_smoke 是限量真实读段的工程联调，不得断言实验合格、科学改善或正式下游可用。
读段数量依次为两个目标重复、两个 IgG 重复。目标峰数和 target_frip 按目标重复顺序排列。
target_frip 是目标片段落在自身候选峰的比例；igg_overlap 顺序为 IgG1 在目标峰1/2中的比例、IgG2 在目标峰1/2中的比例。
replicate_peak_jaccard 描述两个目标重复候选峰的碱基重叠。候选峰数量和 FRiP 本身不能证明改善。
有重跑机会时可以选择 compare_igg，把两份 IgG 作为共同背景运行一次窄峰找峰。第二轮只做复核。
只输出 JSON，字段 action、reason、missing_information。action 为 keep、compare_igg 或 request_evidence。
reason 用中文具体解释。missing_information 是数组，只接受 igv_review 或 quantitative_spikein_design，
且只在 request_evidence 时填写，单独缺 quantitative_spikein_design 不影响本次任务。
不得返回命令、阈值、路径，或杜撰本组没有提供的实验事实。"""


def bounded(value, start, end):
    return type(value) in (int, float) and math.isfinite(value) and start <= value <= end


def validate_payload(payload):
    expected = {'data_source', 'input_scope', 'read_pairs_per_sample', 'round_index', 'experiment_context',
                'baseline', 'with_igg'}
    if (not isinstance(payload, dict) or set(payload) != expected or payload['data_source'] != 'public_h9_summary'
            or payload['experiment_context'] != CONTEXT or payload['input_scope'] not in ('full', 'prefix_smoke')
            or type(payload['round_index']) is not int or payload['round_index'] not in (0, 1)):
        raise ValueError('H9 模型摘要范围无效')
    pairs = payload['read_pairs_per_sample']
    if not isinstance(pairs, list) or len(pairs) != 4 or any(type(x) is not int or not 0 < x <= 20_000_000 for x in pairs):
        raise ValueError('H9 模型摘要读段数无效')
    validate_metrics(payload['baseline'], 0)
    if payload['round_index'] == 1:
        validate_metrics(payload['with_igg'], 2)
    elif payload['with_igg'] is not None:
        raise ValueError('H9 第一轮摘要不能包含重跑结果')


def validate_metrics(metrics, control_count):
    fields = {'peak_mode', 'qvalue', 'control_count', 'peak_counts', 'target_frip', 'igg_overlap', 'replicate_peak_jaccard'}
    if (not isinstance(metrics, dict) or set(metrics) != fields or metrics['peak_mode'] != 'narrow'
            or metrics['qvalue'] != 0.01 or type(metrics['control_count']) is not int
            or metrics['control_count'] != control_count):
        raise ValueError('H9 找峰摘要无效')
    for field, length in (('peak_counts', 2), ('target_frip', 2), ('igg_overlap', 4)):
        values = metrics[field]
        if not isinstance(values, list) or len(values) != length:
            raise ValueError('H9 找峰摘要缺少重复或对照数据')
        if any(not bounded(x, 0, 1 if field != 'peak_counts' else 100_000_000) for x in values):
            raise ValueError('H9 找峰摘要数值无效')
        if field == 'peak_counts' and not all(type(x) is int for x in values):
            raise ValueError('H9 峰计数必须是整数')
    overlap = metrics['replicate_peak_jaccard']
    if overlap is not None and not bounded(overlap, 0, 1):
        raise ValueError('H9 重复峰重叠比例无效')


def validate_decision(value, round_index):
    if (not isinstance(value, dict) or set(value) != {'action', 'reason', 'missing_information'}
            or value['action'] not in ('keep', 'compare_igg', 'request_evidence')
            or not isinstance(value['reason'], str) or not 1 <= len(value['reason'].strip()) <= 1500):
        raise ValueError('H9 模型回复格式无效，未执行调整')
    missing = value['missing_information']
    if (not isinstance(missing, list) or len(missing) > 2
            or any(item not in ('igv_review', 'quantitative_spikein_design') for item in missing)
            or len(set(missing)) != len(missing)
            or ('igv_review' not in missing and value['action'] == 'request_evidence')
            or (missing and value['action'] != 'request_evidence')):
        raise ValueError('H9 模型请求了当前任务不需要的信息')
    if value['action'] == 'compare_igg' and round_index != 0:
        raise ValueError('H9 模型已用完一次重跑机会')
    return value
