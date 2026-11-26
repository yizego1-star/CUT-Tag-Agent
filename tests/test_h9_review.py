import copy
import unittest

from cuttag_web.h9_review import CONTEXT, validate_decision, validate_payload


class H9ReviewTests(unittest.TestCase):
    def test_public_summary_has_fixed_context_and_no_arbitrary_fields(self):
        metrics = {'peak_mode': 'narrow', 'qvalue': 0.01, 'control_count': 0, 'peak_counts': [0, 2],
                   'target_frip': [0, 0.01], 'igg_overlap': [0, 0.01, 0, 0.02], 'replicate_peak_jaccard': 0}
        payload = {'data_source': 'public_h9_summary', 'input_scope': 'prefix_smoke',
                   'read_pairs_per_sample': [100_000] * 4, 'round_index': 0,
                   'experiment_context': copy.deepcopy(CONTEXT), 'baseline': metrics, 'with_igg': None}
        validate_payload(payload)
        payload['experiment_context']['target'] = 'unknown'
        with self.assertRaises(ValueError):
            validate_payload(payload)
        payload['experiment_context'] = CONTEXT
        payload['baseline']['path'] = 'local.bam'
        with self.assertRaises(ValueError):
            validate_payload(payload)

    def test_igg_can_be_compared_once_and_missing_spikein_cannot_block_it(self):
        decision = {'action': 'compare_igg', 'reason': '比较对照片段落峰比例', 'missing_information': []}
        validate_decision(decision, 0)
        with self.assertRaises(ValueError):
            validate_decision(decision, 1)
        decision.update(action='request_evidence', missing_information=['quantitative_spikein_design'])
        with self.assertRaises(ValueError):
            validate_decision(decision, 1)
        decision['missing_information'] = ['igv_review']
        validate_decision(decision, 1)
