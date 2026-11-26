"""Downstream compatibility check with generated sequences, coordinates and annotations only."""
import importlib.util
import json
from pathlib import Path
import pyBigWig
import pysam

spec = importlib.util.spec_from_file_location('worker', '/opt/pipeline-worker.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
out = Path('/output')
(out/'peaks').mkdir(exist_ok=True)
support=out/'support'
(support/'annotation').mkdir(parents=True,exist_ok=True)
with (support/'annotation/gencode.v44.genes.bam-compatible.bed').open('w') as stream:
    for i in range(62700):
        start=10000+i*25
        stream.write('chr1\t%d\t%d\tSYNTHETIC_GENE_%d\t0\t%s\n' % (start,start+500,i,'+' if i%2 else '-'))
with pyBigWig.open(str(out/'coverage.cpm.bw'),'w') as bw:
    bw.addHeader([('chr1',2000000)])
    bw.addEntries(['chr1']*20000,list(range(0,2000000,100)),ends=list(range(100,2000001,100)),values=[float(i%11) for i in range(20000)])
# A synthetic paired-end pileup. Not intended to validate biological inference.
header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000000}]}
records=[]
for i in range(200):
    start=10000+i*2
    for end in (0,1):
        r=pysam.AlignedSegment()
        r.query_name='SYNTHETIC_%d'%i;r.query_sequence='A'*50;r.query_qualities=pysam.qualitystring_to_array('I'*50)
        r.flag=99 if not end else 147;r.reference_id=0;r.reference_start=start+end*100
        r.mapping_quality=40;r.cigarstring='50M';r.next_reference_id=0;r.next_reference_start=start+(1-end)*100
        r.template_length=150 if not end else -150;records.append(r)
with pysam.AlignmentFile(out/'synthetic.bam','wb',header=header) as stream:
    for r in sorted(records,key=lambda a:a.reference_start):stream.write(r)
peak_sets=w.peaks(str(out/'synthetic.bam'),{'chr1':2000000})
info=w.tss()
assert info['genes']==62700 and info['columns']==120
print(json.dumps({'status':'passed','data':'synthetic BAM / BW / annotation only','macs3_modes':list(peak_sets),'tss':info}))
