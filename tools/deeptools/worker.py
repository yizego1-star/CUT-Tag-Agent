"""Container-side fixed coverage step. No network or arbitrary command input."""
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pyBigWig
import pysam


REGION = ("chr1", 1000000, 2000000)
BIN_SIZE = 25


def write_json(name, value):
    Path(name).write_text(json.dumps(value, indent=2) + "\n")


def main():
    shutil.copyfile("/opt/requirements.lock", "/output/requirements.lock")
    if len(sys.argv) not in (2, 3) or Path(sys.argv[1]).name != sys.argv[1]:
        raise ValueError("Expected BAM basename and optional allowed bin size")
    bin_size = int(sys.argv[2]) if len(sys.argv) == 3 else BIN_SIZE
    if bin_size not in (25, 50):
        raise ValueError("Unsupported bin size")
    bam = Path("/input") / sys.argv[1]
    chrom, start, end = REGION
    stat_commands = [["samtools", "quickcheck", "-v", str(bam)], ["samtools", "idxstats", str(bam)]]
    outputs = []
    for command, name in zip(stat_commands, ["samtools-quickcheck.txt", "samtools-idxstats.tsv"]):
        with Path(name).open("w") as stream:
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT)
            invocation = {"command": command, "pid": process.pid, "artifact": name}
            outputs.append(invocation)
            write_json("stat-invocations.json", outputs)
            code = process.wait()
            invocation["exit_code"] = code
            write_json("stat-invocations.json", outputs)
        if code:
            raise ValueError("SAMtools input check/statistics failed; see local log")
    rows = [line.split("\t") for line in Path("samtools-idxstats.tsv").read_text().splitlines()]
    counts = {"mapped_segments": sum(int(row[2]) for row in rows), "unmapped_segments": sum(int(row[3]) for row in rows),
              "reference_count": sum(row[0] != "*" for row in rows), "source_stage": "filtered_deduplicated_bam",
              "quickcheck": "passed", "commands": outputs, "samtools_version": Path("/opt/samtools-version.txt").read_text(errors="replace").splitlines()[0]}
    write_json("bam-stats.json", counts)
    with pysam.AlignmentFile(str(bam), "rb") as source:
        if not source.has_index() or source.get_reference_length(chrom) < end:
            raise ValueError("Required region or BAM index unavailable")
    command = ["bamCoverage", "--bam", str(bam), "--outFileName", "/output/coverage.bw",
               "--outFileFormat", "bigwig", "--binSize", str(bin_size),
               "--region", f"{chrom}:{start}:{end}", "--numberOfProcessors", "2",
               "--normalizeUsing", "None", "--scaleFactor", "1", "--samFlagInclude", "64"]
    before = time.monotonic()
    process = subprocess.Popen(command)
    write_json("tool-start.json", {"command": command, "pid": process.pid,
                                  "tool_version": importlib.metadata.version("deeptools")})
    code = process.wait()
    write_json("tool-exit.json", {"exit_code": code, "elapsed_seconds": round(time.monotonic() - before, 3)})
    if code:
        return code
    with pyBigWig.open("/output/coverage.bw") as track:
        if not track.isBigWig() or track.chroms().get(chrom, 0) < end:
            raise ValueError("Output is not a valid bigWig for the requested region")
        intervals = track.intervals(chrom, start, end) or []
        if not intervals:
            raise ValueError("Output has no intervals in the requested region")
        if any(a < start or b > end or not (float('-inf') < value < float('inf'))
               for a, b, value in intervals):
            raise ValueError("Output intervals are invalid or outside the allowed region")
        values = track.stats(chrom, start, end, nBins=200, type="mean", exact=True)
        result = {"scope": "regional_tool_test", "chrom": chrom, "start": start, "end": end,
                  "bin_size": bin_size, "normalization": "None", "scale_factor": 1,
                  "sam_flag_include": 64, "interval_count": len(intervals),
                  "represented_bases": sum(b - a for a, b, _ in intervals),
                  "nonzero_intervals": sum(value > 0 for _, _, value in intervals),
                  "preview": {"aggregation": "mean", "points": values},
                  "tool_version": importlib.metadata.version("deeptools"),
                  "sha256": hashlib.sha256(Path("coverage.bw").read_bytes()).hexdigest(),
                  "bytes": os.stat("coverage.bw").st_size,
                  "validation": "passed"}
    write_json("result.json", result)
    print("Validated regional bigWig; this is a technical test, not a biological quality verdict.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
