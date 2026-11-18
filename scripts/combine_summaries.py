#!/usr/bin/env python3
"""Combine per-sample EBV summaries and add control-aware classifications."""

import argparse
import csv
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--positive-min-fragments", type=int, default=5)
    parser.add_argument("--positive-min-loci", type=int, default=3)
    parser.add_argument("--background-multiplier", type=float, default=3.0)
    parser.add_argument("--minimum-usable-fragments", type=float, default=100000)
    parser.add_argument("--positive-control-min-fragments", type=int, default=20)
    return parser.parse_args()


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            file_rows = list(reader)
            if len(file_rows) != 1:
                raise ValueError(f"Expected exactly one summary row in {path}")
            rows.append(file_rows[0])
    return rows


def number(row, key, kind=float):
    return kind(row[key])


def rpm_background_ratio(rpm, background_rpm):
    if background_rpm > 0:
        return f"{rpm / background_rpm:.6f}"
    if rpm > 0:
        return "Inf"
    return "NA"


def main():
    args = parse_args()
    rows = read_rows(args.inputs)

    negatives = [row for row in rows if row["control_type"] == "negative"]
    positives = [row for row in rows if row["control_type"] == "positive"]
    negative_control_status = "provided" if negatives else "not_provided"

    background_count = max(
        (number(row, "high_conf_ebv_nonrepeat_fragments", int) for row in negatives),
        default=0,
    )
    background_rpm = max(
        (number(row, "ebv_fragments_per_million") for row in negatives),
        default=0.0,
    )

    if positives:
        positive_control_pass = any(
            number(row, "high_conf_ebv_nonrepeat_fragments", int)
            >= args.positive_control_min_fragments
            and number(row, "distinct_ebv_1kb_bins", int) >= args.positive_min_loci
            for row in positives
        )
        positive_control_status = "pass" if positive_control_pass else "fail"
    else:
        positive_control_pass = None
        positive_control_status = "not_provided"

    for row in rows:
        fragment_count = number(row, "high_conf_ebv_nonrepeat_fragments", int)
        loci = number(row, "distinct_ebv_1kb_bins", int)
        rpm = number(row, "ebv_fragments_per_million")
        usable = number(row, "usable_nonduplicate_fragments_estimate")
        control_type = row["control_type"]

        if control_type == "negative":
            classification = "negative_control"
            reason = "used_to_estimate_background"
        elif control_type == "positive":
            passed = fragment_count >= args.positive_control_min_fragments and loci >= args.positive_min_loci
            classification = "positive_control_pass" if passed else "positive_control_fail"
            reason = "positive_control_acceptance_check"
        elif usable < args.minimum_usable_fragments:
            classification = "qc_failure"
            reason = "insufficient_usable_fragments"
        elif positive_control_pass is False:
            classification = "invalid"
            reason = "positive_control_failed"
        elif fragment_count == 0:
            classification = "negative"
            reason = "no_high_confidence_nonrepeat_ebv_fragments"
        else:
            count_above_background = fragment_count > background_count
            if background_rpm == 0:
                rpm_above_background = True
            else:
                rpm_above_background = rpm >= args.background_multiplier * background_rpm
            if (
                fragment_count >= args.positive_min_fragments
                and loci >= args.positive_min_loci
                and count_above_background
                and rpm_above_background
            ):
                classification = "positive"
                reason = "fragment_count_locus_spread_and_background_criteria_met"
            else:
                classification = "indeterminate"
                reason = "some_ebv_evidence_but_positive_criteria_not_met"

        if (
            control_type == "test"
            and classification in {"negative", "positive", "indeterminate"}
            and (not negatives or not positives)
        ):
            classification = "provisional_" + classification
            missing = []
            if not negatives:
                missing.append("negative")
            if not positives:
                missing.append("positive")
            reason += ";missing_" + "_and_".join(missing) + "_control"

        row["negative_control_background_max_fragments"] = str(background_count)
        row["negative_control_background_max_rpm"] = f"{background_rpm:.6f}"
        row["ebv_rpm_over_negative_background"] = (
            rpm_background_ratio(rpm, background_rpm) if negatives else "NA"
        )
        row["negative_control_status"] = negative_control_status
        row["positive_control_status"] = positive_control_status
        row["evidence_classification"] = classification
        row["classification_reason"] = reason

    fieldnames = list(rows[0].keys())
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
