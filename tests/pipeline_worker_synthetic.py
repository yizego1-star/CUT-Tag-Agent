"""Opt-in container fixture. Only synthetic pairs; no private data or external services."""
import importlib.util
import json
from pathlib import Path
import pysam

spec = importlib.util.spec_from_file_location('worker', '/opt/pipeline-worker.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
out = Path('/output')
(out / 'work').mkdir()
header = {'HD': {'VN': '1.6', 'SO': 'queryname'}, 'SQ': [{'SN': 'chr1', 'LN': 10000}, {'SN': 'chrM', 'LN': 10000}]}


def pair(name, start=100, span=150, mapq=(40, 40), chrom=0, duplicate=False):
    rows = []
    for i in (0, 1):
        r = pysam.AlignedSegment()
        r.query_name, r.query_sequence = name, 'A' * 50
        r.flag = (99 if i == 0 else 147) + (1024 if duplicate else 0)
        r.reference_id, r.reference_start = chrom, start if i == 0 else start + span - 50
        r.mapping_quality = mapq[i]
        r.cigarstring = '50M'
        r.next_reference_id, r.next_reference_start = chrom, start + span - 50 if i == 0 else start
        r.template_length = span if i == 0 else -span
        r.query_qualities = pysam.qualitystring_to_array('I' * 50)
        rows.append(r)
    return rows

with pysam.AlignmentFile(out / 'work/name.bam', 'wb', header=header) as stream:
    for rows in [pair('a-valid'), pair('b-low-mate', mapq=(40, 10)), pair('c-mito', chrom=1),
                 pair('d-long', span=1500), pair('e-orphan')[:1], pair('f-retain-duplicate', start=500, duplicate=True)]:
        for r in rows:
            stream.write(r)
w.filter_pairs()
with pysam.AlignmentFile(out / 'work/filtered.bam', 'rb') as stream:
    kept = [r.query_name for r in stream]
assert kept == ['a-valid', 'a-valid', 'f-retain-duplicate', 'f-retain-duplicate'], kept
peak_sets = {'narrow': {'regions': {'chr1': {'starts': [200], 'ends': [220]}}},
             'broad': {'regions': {'chr1': {'starts': [120, 450], 'ends': [150, 550]}}}}
metrics = w.fragment_metrics('/output/work/filtered.bam', peak_sets)
assert metrics['fragments'] == 2
assert metrics['marked_duplicate_fragments'] == 1
assert metrics['frip']['narrow']['fraction'] == .5
assert metrics['frip']['broad']['fraction'] == 1
assert metrics['fragment_length_counts'][150] == 2
print(json.dumps({'status':'passed','data':'synthetic only','both_mate_mapq':'checked','paired_only':'checked',
                  'chrM_and_length_filter':'checked','duplicate_retention':'checked','fragment_frip':'checked'}))
