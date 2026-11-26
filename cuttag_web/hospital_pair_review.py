"""Fixed, identifier-free model inputs for the registered target and IgG pair."""
from math import isfinite


SCHEMA = 'cuttag_paired_qc_v1'
AUTH_GROUPS = {'settings', 'baseline', 'spikein', 'experiment_context', 'comparison'}
PENDING = {'target_antibody', 'biological_replicates', 'control_design', 'igv_review', 'spikein_quantitation'}
QVALUES = (0.005, 0.02)
SYSTEM = """你负责复核一组 CUT&Tag 目标与匹配 IgG。用中文具体说明已知结果支持什么，以及下一步需要什么。
settings 描述已经完成的候选窄峰试算。baseline 是无 IgG 找峰。target 和 igg 分别代表目标和对照的双端片段数及落在候选峰内的比例。
若选择 compare_igg，本地程序会读取已算好的 IgG 对照结果，并将它作为 comparison 给你复评。两套峰的范围不同，单凭峰数或落峰比例不能宣布哪套科学上更好。结合 IgG 背景、目标信号保留和待做的人工 IGV 检查。匹配关系已知，实验对照设计细节可能待确认。
当前 CPM 用于浏览。Spike-in 是历史过滤计数，标称用量相同也不能消除投入量和计数口径的不确定性。不要将 CPM 与历史 Spike-in 差异归因于标准化，也不要因定量未确认而阻止观察性的候选峰复核。不要编造实验信息。
只输出 JSON，恰好六项：action（keep、request_evidence 或 compare_igg）、bin_size（null）、assessment（中文评估）、reason（中文理由）、missing_information（只填 experiment_context 中确有必要且为 unconfirmed 的字段名）、evidence_used（实际引用的输入组名）。keep 表示保留供探索。request_evidence 需指出至少一个影响下一步的待确认项。只有首轮可选 compare_igg，复评必须引用 comparison。不要输出命令、路径或额外参数。"""


def review_prompt(payload, adaptive=False):
    if adaptive:
        screen = payload['comparison'].get('igv_screen')
        igv = ('comparison.igv_screen 是启动时已准备好的全基因组区域统计。igv_examples_reviewed 是另有人工查看记录的窗口数，0 表示仅做自动统计。'
               'removed 是原无 IgG 峰与任一带 IgG 峰完全不重叠的区间，overlaps_controlled 是其余原峰。'
               'target_higher、igg_higher 分别是对应峰内均值至少为另一轨道 2 倍，similar 表示两者均未达到 2 倍。'
               'target_local_peak 是 target_higher 中目标还至少为等长相邻侧翼平均值 2 倍的个数。'
               '这些是用于阈值探索的 CPM 分组。抗体、重复及实验记录的未知项留到生信验收，勿要求用户补资料后才能继续动作。'
               '如果当前统计没有支持重跑的收益，选择 keep 并解释数值依据。' if screen else '')
        evidence = sorted(AUTH_GROUPS | ({'igv_screen'} if igv else set()))
        pending_action = '' if screen else '、request_evidence'
        phase = ('已经重跑一次。comparison.rerun 是实际新结果。'
                 f'action 只可选 keep（保留原结果）、prefer_rerun（新版在明确指标上有收益）{pending_action}。'
                 'qvalue 填 null。对照写清目标与 IgG 的落峰比例怎样变化和代价。'
                 '如果有 igv_changes，removed_now_overlapping 是原先完全被筛掉后来与新峰重叠的区间；'
                 'overlapping_now_missing 是原先和带 IgG 峰重叠、重跑后完全不重叠的区间。结合这两组的目标局部峰与 IgG 偏强数说明取舍。'
                 'median_igg_cpm_in_rerun_overlap 与 igv_screen.overlaps_controlled.median_igg_cpm 按同一方法统计重跑前后的原峰重叠队列。'
                 if 'rerun' in payload['comparison'] else
                 f'action 可选 keep、compare_igg、tune_igg_qvalue{pending_action}。'
                 'keep 是完成复核后保留当前结果。compare_igg 是比较已有的无 IgG 和带 IgG 候选峰，引用 baseline 与 comparison，解释两套结果的差异。'
                 '只有存在可检验的问题和两套峰的具体数据依据时才选 tune_igg_qvalue。'
                 '调参时 qvalue 只可选 0.005（检验能否进一步减少 IgG 背景）或 0.02（检验能否保留更多目标信号），其他动作填 null。'
                 '理由写出依据、预期变化和要检查的取舍。信息不足或现有结果已足够时直接停止。')
        return (f'你复核已有 CUT&Tag 候选窄峰，两份 BAM 和匹配 IgG 固定。baseline 为无 IgG，comparison 为带 IgG 的原结果，q=0.01。'
                f'只允许改变带 IgG 找峰的 MACS3 q 值，最多真跑一次。{igv}{phase}'
                '不得把峰变多、变少或 FRiP 单项变化当作科学验收。对照目标信号保留与 IgG 背景的取舍。'
                'Spike-in 定量待确认，与本次候选峰阈值探索无关。禁止推断实验事实、输出命令或路径。'
                '只输出一个 JSON 对象，恰好七项：action、bin_size（null）、qvalue、assessment、reason、missing_information、evidence_used。'
                f'assessment 和 reason 都写中文。evidence_used 必须引用 comparison，可附加 igv_screen（如有）。'
                f'组名只可选 {", ".join(evidence)}。'
                f'missing_information 只填这些留待验收的英文原字段名：{", ".join(sorted(pending_fields(payload)))}。无项目时用 []。')
    actions = ['keep', 'request_evidence']
    if 'comparison' not in payload:
        actions.append('compare_igg')
    evidence = sorted(AUTH_GROUPS.intersection(payload))
    return (SYSTEM + f'\n本轮 action 可选 {", ".join(actions)}。evidence_used 只能从 {", ".join(evidence)} 选取组名。'
            f'missing_information 只填这些字段名：{", ".join(sorted(PENDING))}。两个数组的字段名都用原样英文，不填解释。')


def peak_metrics(pair, key):
    return {'peak_count': pair['peak_counts'][key],
            **{role: {'fragments': pair['fragments'][role],
                      'frip': pair['fragment_overlap'][role][key]['fraction']} for role in ('target', 'igg')}}


def make_payload(pair):
    spike = pair['spikein']
    settings = pair['peak_settings']
    if settings['mode'] != 'exploratory_narrow' or settings['blacklist_filter'] != 'not_applied':
        raise ValueError('配对找峰设置与已登记的复核基线不符')
    return {
        'data_source': SCHEMA,
        'settings': {'bin_bp': pair['bin_size'], 'qvalue': settings['qvalue'],
                     'peak_mode': 'narrow', 'normalization': 'CPM', 'blacklist_applied': False},
        'baseline': peak_metrics(pair, 'baseline'),
        'spikein': {'addition_stage': spike['addition_stage'], 'dose_pg': spike['dose_pg'],
                    'cells_per_dose': spike['cells_per_dose'], 'actual_counts_known': spike['actual_cell_counts'] is not None,
                    **{role: spike['historical_fragments'][role] for role in ('target', 'igg')}},
        'experiment_context': {'matched_igg': 'confirmed', **{key: 'unconfirmed' for key in sorted(PENDING)}},
    }


def compare_pair(pair, decision):
    if decision['action'] != 'compare_igg':
        return None
    return {'action': 'compare_igg', **peak_metrics(pair, 'with-igg')}


def adaptive_payload(pair):
    payload = {**make_payload(pair), 'comparison': compare_pair(pair, {'action': 'compare_igg'})}
    if pair.get('igv_screen'):
        payload['comparison']['igv_screen'] = pair['igv_screen']
        payload['experiment_context']['igv_review'] = screen_status(pair['igv_screen'])
    return payload


def pending_fields(payload):
    return {key for key in PENDING if payload['experiment_context'].get(key) == 'unconfirmed'}


def screen_status(screen):
    return 'exploratory_reviewed' if screen['igv_examples_reviewed'] else 'automated_screened'


def validate_igv(screen, peak_count):
    fields = {'ratio_threshold', 'igv_examples_reviewed', 'removed', 'overlaps_controlled'}
    categories = {'target_higher', 'igg_higher', 'similar', 'unscored'}
    if (not isinstance(screen, dict) or set(screen) != fields or screen['ratio_threshold'] != 2
            or type(screen['igv_examples_reviewed']) is not int or not 0 <= screen['igv_examples_reviewed'] <= 100):
        raise ValueError('IGV 探索复核汇总格式无效')
    for group in ('removed', 'overlaps_controlled'):
        row = screen[group]
        counts = categories | {'total', 'target_local_peak'}
        medians = {'median_target_cpm', 'median_igg_cpm'}
        if (not isinstance(row, dict) or set(row) != counts | medians
                or any(type(row[key]) is not int or not 0 <= row[key] <= peak_count for key in counts)
                or any((row[key] is not None if row['total'] == row['unscored'] else
                        type(row[key]) not in (int, float) or not 0 <= row[key] <= 1e6) for key in medians)
                or sum(row[key] for key in categories) != row['total'] or row['target_local_peak'] > row['target_higher']):
            raise ValueError('IGV 峰区间计数不一致')
    if screen['removed']['total'] + screen['overlaps_controlled']['total'] != peak_count:
        raise ValueError('IGV 复核汇总与原始候选峰数不一致')


def validate_igv_changes(changes, screen):
    groups = {'removed_now_overlapping': 'removed', 'overlapping_now_missing': 'overlaps_controlled'}
    categories = {'target_higher', 'igg_higher', 'similar', 'unscored'}
    fields = categories | {'total', 'target_local_peak'}
    median_key = 'median_igg_cpm_in_rerun_overlap'
    if not isinstance(changes, dict) or set(changes) != set(groups) | {median_key}:
        raise ValueError('新阈值的 IGV 区间比较缺失')
    for key, group in groups.items():
        row = changes[key]
        if (not isinstance(row, dict) or set(row) != fields
                or any(type(row[field]) is not int or not 0 <= row[field] <= screen[group][field] for field in fields)
                or row['total'] != sum(row[field] for field in categories) or row['target_local_peak'] > row['target_higher']):
            raise ValueError('新阈值的 IGV 区间统计不一致')
    def scored(row):
        return row['total'] - row['unscored']

    count = scored(screen['overlaps_controlled']) + scored(changes['removed_now_overlapping']) - scored(changes['overlapping_now_missing'])
    value = changes[median_key]
    if (count == 0 and value is not None) or (count > 0 and (type(value) not in (int, float) or not 0 <= value <= 1e6)):
        raise ValueError('新阈值的 IgG 区间信号中位数无效')


def validate_payload(payload, adaptive=False):
    fields = {'data_source', *(AUTH_GROUPS - {'comparison'})}
    if (not isinstance(payload, dict) or set(payload) not in (fields, fields | {'comparison'})
            or payload['data_source'] != SCHEMA):
        raise ValueError('配对复核摘要字段不符合约定')
    if payload['settings'] != {'bin_bp': 25, 'qvalue': 0.01, 'peak_mode': 'narrow', 'normalization': 'CPM', 'blacklist_applied': False}:
        raise ValueError('配对复核参数与已校验试算不符')
    context = {'matched_igg': 'confirmed', **dict.fromkeys(PENDING, 'unconfirmed')}
    metrics = payload['baseline']
    validate_metrics(metrics)
    if adaptive and 'comparison' not in payload:
        raise ValueError('自动调参需要已完成的带 IgG 对照基线')
    if 'comparison' in payload:
        comparison = payload['comparison']
        if not isinstance(comparison, dict) or comparison.get('action') != 'compare_igg':
            raise ValueError('配对复评类型不符合约定')
        excluded = {'action', 'rerun', 'igv_screen'} if adaptive else {'action'}
        validate_metrics({key: value for key, value in comparison.items() if key not in excluded})
        if adaptive and 'igv_screen' in comparison:
            validate_igv(comparison['igv_screen'], metrics['peak_count'])
            context['igv_review'] = screen_status(comparison['igv_screen'])
        if any(comparison[role]['fragments'] != metrics[role]['fragments'] for role in ('target', 'igg')):
            raise ValueError('两套候选峰对应的片段总数不一致')
        if adaptive and 'rerun' in comparison:
            rerun = comparison['rerun']
            if not isinstance(rerun, dict) or rerun.get('qvalue') not in QVALUES:
                raise ValueError('新候选峰的调参阈值无效')
            excluded_rerun = {'qvalue'}
            if 'igv_screen' in comparison:
                validate_igv_changes(rerun.get('igv_changes'), comparison['igv_screen'])
                excluded_rerun.add('igv_changes')
            validate_metrics({key: value for key, value in rerun.items() if key not in excluded_rerun})
            if any(rerun[role]['fragments'] != metrics[role]['fragments'] for role in ('target', 'igg')):
                raise ValueError('新候选峰与原峰使用的配对片段总数不一致')
    if payload['experiment_context'] != context:
        raise ValueError('配对复核实验状态不符合约定')
    spike = payload['spikein']
    if (not isinstance(spike, dict) or set(spike) != {'addition_stage', 'dose_pg', 'cells_per_dose', 'actual_counts_known', 'target', 'igg'}
            or spike['addition_stage'] != 'after_fragmentation' or spike['dose_pg'] != 1 or spike['cells_per_dose'] != 200000
            or type(spike['actual_counts_known']) is not bool
            or any(type(spike[role]) is not int or not 0 <= spike[role] <= 10**9 for role in ('target', 'igg'))):
        raise ValueError('配对复核 Spike-in 摘要不符合约定')


def validate_metrics(metrics):
    if (not isinstance(metrics, dict) or set(metrics) != {'peak_count', 'target', 'igg'}
            or type(metrics['peak_count']) is not int or not 0 <= metrics['peak_count'] <= 10**9):
        raise ValueError('配对复核的候选峰计数无效')
    for role in ('target', 'igg'):
        values = metrics[role]
        if (not isinstance(values, dict) or set(values) != {'fragments', 'frip'}
                or type(values['fragments']) is not int or not 0 < values['fragments'] <= 10**9
                or type(values['frip']) not in (int, float) or not 0 <= values['frip'] <= 1 or not isfinite(values['frip'])):
            raise ValueError('配对复核的片段计数或落峰比例无效')


def resolve_igv_evidence(keys, payload):
    """Accept references to existing subfields and store their stable display group."""
    resolved = []
    for key in keys:
        if not isinstance(key, str):
            raise ValueError('配对复核只允许引用已提供的证据')
        if key == 'igv_screen':
            key = 'comparison.igv_screen'
        path = key.split('.')
        value = payload
        if path[0] not in AUTH_GROUPS:
            raise ValueError('配对复核只允许引用已提供的证据')
        for part in path:
            if not isinstance(value, dict) or part not in value:
                raise ValueError('配对复核只允许引用已提供的证据')
            value = value[part]
        group = 'igv_screen' if path[:2] == ['comparison', 'igv_screen'] else path[0]
        if group not in resolved:
            resolved.append(group)
    return resolved


def validate_decision(value, payload, adaptive=False):
    fields = {'action', 'bin_size', 'assessment', 'reason', 'missing_information', 'evidence_used'}
    rerun = payload.get('comparison', {}).get('rerun') if adaptive else None
    actions = {'keep', 'request_evidence'}
    if adaptive:
        fields.add('qvalue')
        actions.add('prefer_rerun' if rerun else 'tune_igg_qvalue')
        if not rerun:
            actions.add('compare_igg')
        if 'igv_screen' in payload.get('comparison', {}):
            actions.discard('request_evidence')
    else:
        actions.add('compare_igg')
    if (not isinstance(value, dict) or set(value) != fields or value['bin_size'] is not None
            or not isinstance(value['action'], str) or value['action'] not in actions
            or any(not isinstance(value[key], str) or not 1 <= len(value[key].strip()) <= 2000 for key in ('assessment', 'reason'))):
        raise ValueError('配对复核的模型回复格式无效')
    missing, evidence = value['missing_information'], value['evidence_used']
    pending = pending_fields(payload)
    if (not isinstance(missing, list) or any(not isinstance(key, str) or key not in pending for key in missing)
            or len(set(missing)) != len(missing) or (value['action'] == 'request_evidence' and not missing)):
        raise ValueError('模型所请求的信息超出当前待确认范围')
    available = AUTH_GROUPS.intersection(payload)
    if adaptive and 'igv_screen' in payload.get('comparison', {}) and isinstance(evidence, list):
        evidence = resolve_igv_evidence(evidence, payload)
        available.add('igv_screen')
        value = {**value, 'evidence_used': evidence}
    if (not isinstance(evidence, list) or not evidence or any(not isinstance(key, str) or key not in available for key in evidence)
            or len(set(evidence)) != len(evidence)):
        raise ValueError('配对复核只允许引用已提供的证据')
    cited_comparison = 'comparison' in evidence or ('igv_screen' in evidence and not rerun)
    if ('comparison' in payload and (not cited_comparison or (not adaptive and value['action'] == 'compare_igg'))
            or value['action'] == 'compare_igg' and 'baseline' not in evidence):
        raise ValueError('比较仅可执行一次，且必须引用已提供的峰证据')
    if adaptive:
        tuning = value['action'] == 'tune_igg_qvalue'
        if (tuning and value['qvalue'] not in QVALUES) or (not tuning and value['qvalue'] is not None):
            raise ValueError('调参只能选择允许的 q 值，结束动作不可再改参数')
        if value['action'] == 'prefer_rerun':
            previous = payload['comparison']
            moves_in_goal_direction = (rerun['igg']['frip'] < previous['igg']['frip'] if rerun['qvalue'] < .01 else
                                       rerun['target']['frip'] > previous['target']['frip'])
            if rerun['qvalue'] > .01 and 'igv_changes' in rerun:
                moves_in_goal_direction = moves_in_goal_direction and rerun['igv_changes']['removed_now_overlapping']['target_higher'] > 0
            if not moves_in_goal_direction or not rerun['target']['frip']:
                raise ValueError('实际统计没有支持新阈值的预期收益')
    return value
