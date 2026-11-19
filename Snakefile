import csv
import re
import shlex


configfile: "config/config.yaml"


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def load_samples(path):
    required = {"sample", "r1", "r2"}
    units = {}
    sample_units = {}
    controls = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"samples.tsv is missing columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            sample = row["sample"].strip()
            if not sample or sample.startswith("#"):
                continue
            if not SAFE_ID.fullmatch(sample):
                raise ValueError(
                    f"Invalid sample name {sample!r}; use only letters, numbers, ., _, or -"
                )
            sample_units.setdefault(sample, [])
            unit = row.get("unit", "").strip() or f"unit{len(sample_units[sample]) + 1}"
            if not SAFE_ID.fullmatch(unit):
                raise ValueError(
                    f"Invalid unit name {unit!r}; use only letters, numbers, ., _, or -"
                )
            key = (sample, unit)
            if key in units:
                raise ValueError(f"Duplicate sample/unit combination: {sample}/{unit}")
            r1 = row["r1"].strip()
            r2 = row["r2"].strip()
            if not r1 or not r2:
                raise ValueError(f"Missing R1 or R2 path for {sample}/{unit}")
            units[key] = {"r1": r1, "r2": r2}
            sample_units[sample].append(unit)

            raw_control = row.get("control_type", "").strip().lower()
            control = raw_control or controls.get(sample, "test")
            if control not in {"test", "negative", "positive"}:
                raise ValueError(
                    f"control_type for {sample} must be test, negative, or positive"
                )
            if sample in controls and controls[sample] != control:
                raise ValueError(
                    f"All rows for {sample} must have the same control_type"
                )
            controls[sample] = control
    if not units:
        raise ValueError("No samples were found in samples.tsv")
    return units, sample_units, controls


UNITS, SAMPLE_UNITS, CONTROL_TYPES = load_samples(config["samples"])
SAMPLE_NAMES = list(SAMPLE_UNITS)
BWA_INDEX = config["bwa_index"]
BWA_INDEX_FILES = [
    f"{BWA_INDEX}{suffix}" for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
]

BWA_THREADS = int(config.get("bwa_threads", 6))
SORT_THREADS = int(config.get("sort_threads", 2))
SORT_MEMORY = str(config.get("sort_memory_per_thread", "1G"))
SAMTOOLS_THREADS = int(config.get("samtools_threads", 2))
MARKDUP_BARCODE_OPTION = (
    "--barcode-name" if config.get("barcode_aware_duplicates", True) else ""
)


wildcard_constraints:
    sample=r"[A-Za-z0-9_.-]+",
    unit=r"[A-Za-z0-9_.-]+",


rule all:
    input:
        "results/ebv_summary.tsv",
        expand("results/bam/{sample}.hg38_EBV.bam.bai", sample=SAMPLE_NAMES),
        expand("results/evidence/{sample}.ebv_candidates.tsv", sample=SAMPLE_NAMES),


rule align_bam:
    input:
        r1=lambda wc: UNITS[(wc.sample, wc.unit)]["r1"],
        r2=lambda wc: UNITS[(wc.sample, wc.unit)]["r2"],
        index_files=BWA_INDEX_FILES,
    output:
        temp("results/intermediate/units/{sample}/{unit}.unsorted.bam"),
    log:
        bwa="results/logs/units/{sample}/{unit}.bwa.log",
        bam="results/logs/units/{sample}/{unit}.bam_conversion.log",
    params:
        index=lambda wc, input: str(input.index_files[0])[:-4],
        bwa_threads=BWA_THREADS,
    threads: BWA_THREADS + 1
    conda: "environment.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/intermediate/units/{wildcards.sample} \
          results/logs/units/{wildcards.sample}
        bwa mem -t {params.bwa_threads} {params.index:q} {input.r1:q} {input.r2:q} \
          2> {log.bwa:q} \
        | samtools view -@ 1 -b -o {output:q} - 2> {log.bam:q}
        """


rule name_sort_unit:
    input:
        "results/intermediate/units/{sample}/{unit}.unsorted.bam",
    output:
        temp("results/intermediate/units/{sample}/{unit}.namesort.bam"),
    log:
        "results/logs/units/{sample}/{unit}.namesort.log",
    params:
        sort_memory=SORT_MEMORY,
    threads: SORT_THREADS
    conda: "environment.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/tmp
        sort_tmpdir=$(mktemp -d \
          results/tmp/{wildcards.sample}.{wildcards.unit}.namesort.XXXXXX)
        trap 'rm -rf "$sort_tmpdir"' EXIT
        samtools sort -n -@ {threads} -m {params.sort_memory} \
          -T "$sort_tmpdir/sort" \
          -o {output:q} {input:q} 2> {log:q}
        """


def sample_unit_bams(wildcards):
    return [
        f"results/intermediate/units/{wildcards.sample}/{unit}.namesort.bam"
        for unit in SAMPLE_UNITS[wildcards.sample]
    ]


rule merge_name_sorted_units:
    input:
        sample_unit_bams,
    output:
        temp("results/intermediate/{sample}.namesort.bam"),
    log:
        "results/logs/{sample}.merge_units.log",
    threads: SAMTOOLS_THREADS
    conda: "environment.yml"
    shell:
        "samtools merge -n -f -@ {threads} {output:q} {input:q} 2> {log:q}"


rule fixmate:
    input:
        "results/intermediate/{sample}.namesort.bam",
    output:
        temp("results/intermediate/{sample}.fixmate.bam"),
    log:
        "results/logs/{sample}.fixmate.log",
    threads: SAMTOOLS_THREADS
    conda: "environment.yml"
    shell:
        "samtools fixmate -@ {threads} -m {input:q} {output:q} 2> {log:q}"


rule coordinate_sort:
    input:
        "results/intermediate/{sample}.fixmate.bam",
    output:
        temp("results/intermediate/{sample}.positionsort.bam"),
    log:
        "results/logs/{sample}.positionsort.log",
    params:
        sort_memory=SORT_MEMORY,
    threads: SORT_THREADS
    conda: "environment.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/tmp
        sort_tmpdir=$(mktemp -d \
          results/tmp/{wildcards.sample}.positionsort.XXXXXX)
        trap 'rm -rf "$sort_tmpdir"' EXIT
        samtools sort -@ {threads} -m {params.sort_memory} \
          -T "$sort_tmpdir/sort" \
          -o {output:q} {input:q} 2> {log:q}
        """


rule mark_duplicates:
    input:
        "results/intermediate/{sample}.positionsort.bam",
    output:
        "results/bam/{sample}.hg38_EBV.bam",
    log:
        "results/logs/{sample}.markdup.log",
    params:
        barcode_option=MARKDUP_BARCODE_OPTION,
    threads: SAMTOOLS_THREADS
    conda: "environment.yml"
    shell:
        r"""
        mkdir -p results/bam
        samtools markdup -@ {threads} -s {params.barcode_option} \
          {input:q} {output:q} 2> {log:q}
        """


rule index_bam:
    input:
        "results/bam/{sample}.hg38_EBV.bam",
    output:
        "results/bam/{sample}.hg38_EBV.bam.bai",
    log:
        "results/logs/{sample}.index.log",
    threads: SAMTOOLS_THREADS
    conda: "environment.yml"
    shell:
        "samtools index -@ {threads} {input:q} {output:q} 2> {log:q}"


rule validate_bam:
    input:
        bam="results/bam/{sample}.hg38_EBV.bam",
        bai="results/bam/{sample}.hg38_EBV.bam.bai",
    output:
        touch("results/qc/{sample}.quickcheck.ok"),
    log:
        "results/logs/{sample}.quickcheck.log",
    conda: "environment.yml"
    shell:
        r"""
        mkdir -p results/qc
        samtools quickcheck -v {input.bam:q} > {log:q} 2>&1
        """


rule summarize_sample:
    input:
        bam="results/bam/{sample}.hg38_EBV.bam",
        bai="results/bam/{sample}.hg38_EBV.bam.bai",
        checked="results/qc/{sample}.quickcheck.ok",
    output:
        summary="results/evidence/{sample}.summary.tsv",
        candidates="results/evidence/{sample}.ebv_candidates.tsv",
    log:
        "results/logs/{sample}.summary.log",
    params:
        control_type=lambda wc: CONTROL_TYPES[wc.sample],
        input_fastq_pairs=lambda wc: len(SAMPLE_UNITS[wc.sample]),
        ebv_contig=config["ebv_contig"],
        min_mapq=int(config.get("min_mapq", 30)),
        min_aligned_bases=int(config.get("min_aligned_bases_per_mate", 30)),
        max_nm=int(config.get("max_edit_distance_per_mate", 5)),
        exclude_option=(
            f"--exclude-bed {shlex.quote(str(config['ebv_exclude_bed']))}"
            if config.get("ebv_exclude_bed")
            else ""
        ),
    threads: SAMTOOLS_THREADS
    conda: "environment.yml"
    shell:
        r"""
        mkdir -p results/evidence
        python3 scripts/summarize_ebv.py \
          --sample {wildcards.sample} \
          --control-type {params.control_type} \
          --input-fastq-pairs {params.input_fastq_pairs} \
          --bam {input.bam:q} \
          --ebv-contig {params.ebv_contig:q} \
          --min-mapq {params.min_mapq} \
          --min-aligned-bases {params.min_aligned_bases} \
          --max-nm {params.max_nm} \
          {params.exclude_option} \
          --samtools-threads {threads} \
          --summary {output.summary:q} \
          --candidates {output.candidates:q} \
          > {log:q} 2>&1
        """


rule combine_summaries:
    input:
        expand("results/evidence/{sample}.summary.tsv", sample=SAMPLE_NAMES),
    output:
        "results/ebv_summary.tsv",
    log:
        "results/logs/combine_summaries.log",
    params:
        positive_min_fragments=int(config.get("positive_min_fragments", 5)),
        positive_min_loci=int(config.get("positive_min_distinct_loci", 3)),
        background_multiplier=float(config.get("background_rpm_multiplier", 3.0)),
        minimum_usable_fragments=int(config.get("minimum_usable_fragments", 100000)),
        positive_control_min_fragments=int(config.get("positive_control_min_fragments", 20)),
    conda: "environment.yml"
    shell:
        r"""
        python3 scripts/combine_summaries.py \
          --inputs {input:q} \
          --output {output:q} \
          --positive-min-fragments {params.positive_min_fragments} \
          --positive-min-loci {params.positive_min_loci} \
          --background-multiplier {params.background_multiplier} \
          --minimum-usable-fragments {params.minimum_usable_fragments} \
          --positive-control-min-fragments {params.positive_control_min_fragments} \
          > {log:q} 2>&1
        """
