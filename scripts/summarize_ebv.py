#!/usr/bin/env python3
"""Summarize human and EBV evidence from a coordinate-sorted BAM.

Only the Python standard library and the samtools executable are used.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict


SECONDARY = 0x100
SUPPLEMENTARY = 0x800
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
PROPER_PAIR = 0x2
DUPLICATE = 0x400
QCFAIL = 0x200

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--control-type", default="test")
    parser.add_argument("--input-fastq-pairs", type=int, default=1)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ebv-contig", required=True)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-aligned-bases", type=int, default=30)
    parser.add_argument("--max-nm", type=int, default=5)
    parser.add_argument("--exclude-bed", default="")
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--candidates", required=True)
    return parser.parse_args()


def cigar_metrics(cigar):
    aligned_query = 0
    reference_span = 0
    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in {"M", "I", "=", "X"}:
            aligned_query += length
        if operation in {"M", "D", "N", "=", "X"}:
            reference_span += length
    return aligned_query, reference_span


def parse_tags(fields):
    tags = {}
    for field in fields[11:]:
        pieces = field.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]
    return tags


def load_exclusions(path, contig):
    intervals = []
    if not path:
        return intervals
    if not os.path.exists(path):
        raise FileNotFoundError(f"EBV exclusion BED does not exist: {path}")
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED line {line_number}: {line.rstrip()}")
            if fields[0] == contig:
                intervals.append((int(fields[1]), int(fields[2])))
    intervals.sort()
    return intervals


def overlaps_exclusion(start0, end0, exclusions):
    return any(start0 < interval_end and end0 > interval_start
               for interval_start, interval_end in exclusions)


def union_covered_bases(intervals, contig_length):
    """Count bases covered by the union of 1-based inclusive intervals."""
    clipped = sorted(
        (max(1, start), min(contig_length, end))
        for start, end in intervals
        if start <= contig_length and end >= 1
    )
    if not clipped:
        return 0
    covered = 0
    current_start, current_end = clipped[0]
    for start, end in clipped[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    return covered + current_end - current_start + 1


def record_from_fields(fields):
    flag = int(fields[1])
    position1 = int(fields[3])
    mapq = int(fields[4])
    aligned_query, reference_span = cigar_metrics(fields[5])
    tags = parse_tags(fields)
    nm = int(tags["NM"]) if "NM" in tags else None
    end1 = position1 + max(reference_span, 1) - 1
    return {
        "qname": fields[0],
        "flag": flag,
        "rname": fields[2],
        "position1": position1,
        "end1": end1,
        "mapq": mapq,
        "cigar": fields[5],
        "rnext": fields[6],
        "pnext": int(fields[7]),
        "tlen": int(fields[8]),
        "aligned_query": aligned_query,
        "nm": nm,
        "as": tags.get("AS", ""),
        "duplicate": bool(flag & DUPLICATE),
    }


def mate_passes(record, args):
    if not (record["flag"] & PROPER_PAIR):
        return False, "not_proper_pair"
    if record["flag"] & (UNMAPPED | MATE_UNMAPPED | QCFAIL):
        return False, "unmapped_or_qcfail"
    if record["mapq"] < args.min_mapq:
        return False, "low_mapq"
    if record["aligned_query"] < args.min_aligned_bases:
        return False, "short_alignment"
    if args.max_nm >= 0 and (record["nm"] is None or record["nm"] > args.max_nm):
        return False, "high_or_missing_nm"
    return True, "pass"


def write_tsv(path, fieldnames, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bam_contigs(bam):
    """Return {contig: length} from the BAM header, failing on bad headers."""
    command = ["samtools", "view", "-H", bam]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    contigs = {}
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(
            field.split(":", 1)
            for field in line.split("\t")[1:]
            if ":" in field
        )
        if "SN" in fields and "LN" in fields:
            contigs[fields["SN"]] = int(fields["LN"])
    return contigs


def main():
    args = parse_args()
    exclusions = load_exclusions(args.exclude_bed, args.ebv_contig)
    contigs = bam_contigs(args.bam)
    if args.ebv_contig not in contigs:
        preview = ", ".join(list(contigs)[:8]) or "no @SQ records"
        raise ValueError(
            f"EBV contig {args.ebv_contig!r} is absent from the BAM header; "
            f"first header contigs: {preview}"
        )

    command = ["samtools", "view", "-@", str(args.samtools_threads), args.bam]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)

    metrics = defaultdict(int)
    ebv_records = defaultdict(list)

    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & (SECONDARY | SUPPLEMENTARY):
            continue

        metrics["total_primary_reads"] += 1
        mapped = not bool(flag & UNMAPPED)
        if mapped:
            metrics["mapped_primary_reads"] += 1
            if fields[2] == args.ebv_contig:
                metrics["ebv_mapped_reads"] += 1
                ebv_records[fields[0]].append(record_from_fields(fields))
            else:
                metrics["human_mapped_reads"] += 1
        else:
            metrics["unmapped_primary_reads"] += 1

        usable_record = (
            mapped
            and bool(flag & PROPER_PAIR)
            and not bool(flag & (MATE_UNMAPPED | QCFAIL | DUPLICATE))
            and int(fields[4]) >= args.min_mapq
        )
        if usable_record:
            metrics["usable_nonduplicate_read_records"] += 1

    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"samtools view failed with exit status {return_code}")

    candidate_rows = []
    high_conf_intervals = set()
    high_conf_bins = set()
    high_conf_nonduplicate = 0
    high_conf_nonrepeat = 0
    high_conf_including_duplicates = 0
    proper_ebv_pairs = 0

    for qname in sorted(ebv_records):
        records = ebv_records[qname]
        reasons = []
        if len(records) != 2:
            reasons.append(f"expected_2_primary_ebv_records_found_{len(records)}")
        else:
            if all(record["flag"] & PROPER_PAIR for record in records):
                proper_ebv_pairs += 1
            for record in records:
                passed, reason = mate_passes(record, args)
                if not passed:
                    reasons.append(reason)

        fragment_start1 = min(record["position1"] for record in records)
        fragment_end1 = max(record["end1"] for record in records)
        duplicate = any(record["duplicate"] for record in records)
        excluded = overlaps_exclusion(fragment_start1 - 1, fragment_end1, exclusions)
        pair_pass = len(records) == 2 and not reasons

        if pair_pass:
            high_conf_including_duplicates += 1
            if not duplicate:
                high_conf_nonduplicate += 1
                if not excluded:
                    high_conf_nonrepeat += 1
                    high_conf_intervals.add((fragment_start1, fragment_end1))
                    high_conf_bins.add((fragment_start1 - 1) // 1000)

        if reasons:
            status = "reject:" + ",".join(sorted(set(reasons)))
        elif duplicate:
            status = "reject:duplicate"
        elif excluded:
            status = "reject:excluded_region"
        else:
            status = "retain"

        candidate_rows.append({
            "sample": args.sample,
            "qname": qname,
            "status": status,
            "record_count": len(records),
            "fragment_start_1based": fragment_start1,
            "fragment_end_1based": fragment_end1,
            "mapq": ",".join(str(record["mapq"]) for record in records),
            "cigar": ",".join(record["cigar"] for record in records),
            "aligned_query_bases": ",".join(str(record["aligned_query"]) for record in records),
            "nm": ",".join("" if record["nm"] is None else str(record["nm"]) for record in records),
            "alignment_score": ",".join(str(record["as"]) for record in records),
            "duplicate": str(duplicate).lower(),
            "excluded_region": str(excluded).lower(),
        })

    total_primary_reads = metrics["total_primary_reads"]
    total_read_pairs = total_primary_reads / 2.0
    usable_fragments_estimate = metrics["usable_nonduplicate_read_records"] / 2.0
    rpm = (1_000_000.0 * high_conf_nonrepeat / usable_fragments_estimate
           if usable_fragments_estimate else 0.0)
    ebv_covered_bases = union_covered_bases(
        high_conf_intervals, contigs[args.ebv_contig]
    )
    ebv_genome_breadth = ebv_covered_bases / contigs[args.ebv_contig]

    summary = {
        "sample": args.sample,
        "control_type": args.control_type,
        "input_fastq_pairs": args.input_fastq_pairs,
        "bam": args.bam,
        "ebv_contig": args.ebv_contig,
        "ebv_contig_length": contigs[args.ebv_contig],
        "total_primary_reads": total_primary_reads,
        "total_read_pairs": f"{total_read_pairs:.1f}",
        "mapped_primary_reads": metrics["mapped_primary_reads"],
        "mapped_read_rate": f"{metrics['mapped_primary_reads'] / total_primary_reads:.6f}" if total_primary_reads else "0",
        "human_mapped_reads": metrics["human_mapped_reads"],
        "ebv_mapped_reads": metrics["ebv_mapped_reads"],
        "unmapped_primary_reads": metrics["unmapped_primary_reads"],
        "ebv_mapped_fragments_any": len(ebv_records),
        "ebv_proper_pair_fragments": proper_ebv_pairs,
        "high_conf_ebv_fragments_including_duplicates": high_conf_including_duplicates,
        "high_conf_ebv_nonduplicate_fragments": high_conf_nonduplicate,
        "high_conf_ebv_nonrepeat_fragments": high_conf_nonrepeat,
        "distinct_ebv_fragment_intervals": len(high_conf_intervals),
        "distinct_ebv_1kb_bins": len(high_conf_bins),
        "ebv_covered_bases": ebv_covered_bases,
        "ebv_genome_breadth_fraction": f"{ebv_genome_breadth:.8f}",
        "usable_nonduplicate_read_records": metrics["usable_nonduplicate_read_records"],
        "usable_nonduplicate_fragments_estimate": f"{usable_fragments_estimate:.1f}",
        "ebv_fragments_per_million": f"{rpm:.6f}",
        "min_mapq": args.min_mapq,
        "min_aligned_bases_per_mate": args.min_aligned_bases,
        "max_edit_distance_per_mate": args.max_nm,
        "exclusion_bed_used": args.exclude_bed or "none",
    }

    write_tsv(args.candidates, list(candidate_rows[0].keys()) if candidate_rows else [
        "sample", "qname", "status", "record_count", "fragment_start_1based",
        "fragment_end_1based", "mapq", "cigar", "aligned_query_bases", "nm",
        "alignment_score", "duplicate", "excluded_region"
    ], candidate_rows)
    write_tsv(args.summary, list(summary.keys()), [summary])


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
