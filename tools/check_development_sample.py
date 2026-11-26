#!/usr/bin/env python3
"""Check a locally imported sample. Never connects to or modifies the HPC."""

import argparse
import gzip
import hashlib
import json
import re
import struct
from pathlib import Path


BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def md5sum(path):
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fastq_first_record(path):
    with gzip.open(path, "rb") as stream:
        lines = [stream.readline(65536).rstrip(b"\r\n") for _ in range(4)]
    if not (lines[0].startswith(b"@") and lines[2].startswith(b"+")):
        raise ValueError("First FASTQ record has invalid identifier or separator")
    if not lines[1] or len(lines[1]) != len(lines[3]):
        raise ValueError("First FASTQ sequence and quality lengths do not match")
    identifier = lines[0].split()[0]
    if identifier.endswith((b"/1", b"/2")):
        identifier = identifier[:-2]
    return identifier, len(lines[1])


def bam_header(path):
    with gzip.open(path, "rb") as stream:
        if stream.read(4) != b"BAM\x01":
            raise ValueError("Invalid BAM magic")
        length_bytes = stream.read(4)
        if len(length_bytes) != 4:
            raise ValueError("Truncated BAM header")
        text_length = struct.unpack("<i", length_bytes)[0]
        if not 0 <= text_length <= 10 * 1024 * 1024:
            raise ValueError("Unexpected BAM header size")
        if len(stream.read(text_length)) != text_length:
            raise ValueError("Truncated BAM header text")
        reference_bytes = stream.read(4)
        if len(reference_bytes) != 4:
            raise ValueError("Missing BAM reference count")
        references = struct.unpack("<i", reference_bytes)[0]
        if not 0 < references <= 1000000:
            raise ValueError("Unexpected BAM reference count")
    with path.open("rb") as stream:
        stream.seek(-len(BGZF_EOF), 2)
        eof_ok = stream.read() == BGZF_EOF
    if not eof_ok:
        raise ValueError("BAM is missing the standard BGZF EOF block")
    return references


def verify(manifest_path):
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    root = manifest_path.parent.parent.resolve()
    result = {
        "sample_id": manifest["sample_id"],
        "scope": "local transfer and basic format checks; no scientific QC verdict",
        "checks": [],
        "errors": [],
    }
    files = []
    for item in manifest["files"]:
        path = (root / item["local_relative_path"]).resolve()
        if root not in path.parents:
            result["errors"].append("Manifest path escapes the sample directory")
            continue
        size = path.stat().st_size if path.is_file() else None
        ok = size == item["expected_bytes"]
        files.append((item, path))
        result["checks"].append({
            "check": "file_size", "file": item["local_relative_path"],
            "expected_bytes": item["expected_bytes"], "actual_bytes": size,
            "passed": ok,
        })
        if not ok:
            result["errors"].append("Missing or incomplete: " + item["local_relative_path"])
    if result["errors"]:
        result["status"] = "incomplete"
        return result

    first_records = []
    for item, path in files:
        relative = item["local_relative_path"]
        try:
            if relative.endswith(".fq.gz"):
                sidecar = Path(str(path) + ".md5")
                expected = sidecar.read_text().split()[0].lower()
                if not re.fullmatch(r"[0-9a-f]{32}", expected):
                    raise ValueError("Invalid MD5 sidecar")
                observed = md5sum(path)
                if observed != expected:
                    raise ValueError("MD5 does not match the source sidecar")
                result["checks"].append({"check": "source_md5", "file": relative, "passed": True})
                identifier, read_length = fastq_first_record(path)
                first_records.append(identifier)
                result["checks"].append({
                    "check": "first_fastq_record", "file": relative,
                    "read_length": read_length, "passed": True,
                })
            elif relative.endswith(".bam"):
                references = bam_header(path)
                with Path(str(path) + ".bai").open("rb") as index:
                    if index.read(4) != b"BAI\x01":
                        raise ValueError("Invalid BAI magic")
                    count_bytes = index.read(4)
                    if len(count_bytes) != 4 or struct.unpack("<i", count_bytes)[0] != references:
                        raise ValueError("BAM and BAI reference counts do not match")
                result["checks"].append({
                    "check": "bam_header_eof_and_index_reference_count",
                    "file": relative, "references": references, "passed": True,
                })
            elif relative == "qc/fastp.json":
                data = json.loads(path.read_text(encoding="utf-8"))
                summary = data["summary"]
                for name in ("before_filtering", "after_filtering"):
                    if not isinstance(summary[name]["total_reads"], int):
                        raise ValueError("Invalid fastp total_reads")
                result["checks"].append({"check": "fastp_summary_schema", "passed": True})
        except (OSError, EOFError, ValueError, KeyError, IndexError, struct.error) as exc:
            result["errors"].append(relative + ": " + str(exc))
    if len(first_records) == 2:
        same = first_records[0] == first_records[1]
        result["checks"].append({"check": "first_mate_identifiers_match", "passed": same})
        if not same:
            result["errors"].append("First R1/R2 identifiers do not match")
    result["status"] = "passed" if not result["errors"] else "failed"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.manifest)
    serialized = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    print(json.dumps({
        "sample_id": result["sample_id"], "status": result["status"],
        "check_count": len(result["checks"]), "errors": result["errors"],
    }, ensure_ascii=False))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
