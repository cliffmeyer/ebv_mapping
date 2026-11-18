# Minimal hg38 + EBV ATAC-seq screen

This workflow uses only BWA, SAMtools, Snakemake, and Python's standard library.
It accepts untrimmed paired-end FASTQs, competitively maps them to a prebuilt
hg38+EBV BWA index, marks duplicates, and writes a control-aware TSV summary.

## Inputs

1. Build a combined hg38+EBV FASTA and BWA index once:

   ```bash
   bwa index -p /path/to/hg38_NC_007605.1 combined_hg38_EBV.fa
   ```

2. Edit `config/config.yaml` with the index prefix and exact EBV contig name.

3. Edit `config/samples.tsv`:

   ```text
   sample  unit  r1  r2  control_type
   sample1 L001  /data/sample1_L001_R1.fastq.gz /data/sample1_L001_R2.fastq.gz test
   sample1 L002  /data/sample1_L002_R1.fastq.gz /data/sample1_L002_R2.fastq.gz test
   sample2 L001  /data/sample2_L001_R1.fastq.gz /data/sample2_L001_R2.fastq.gz test
   neg1    L001  /data/neg1_R1.fastq.gz         /data/neg1_R2.fastq.gz         negative
   pos1    L001  /data/pos1_R1.fastq.gz         /data/pos1_R2.fastq.gz         positive
   ```

`control_type` must be `test`, `negative`, or `positive`. Controls are optional,
but classifications without them are provisional. Add one row per paired FASTQ
set. Repeat the sample name for every lane or sequencing unit belonging to that
sample; `unit` must be unique within a sample. If the `unit` column is omitted,
the workflow assigns `unit1`, `unit2`, and so on in row order.

Each FASTQ pair is aligned independently, allowing Snakemake to process pairs
and samples in parallel. Name-sorted unit BAMs are then merged per sample before
mate fixing and duplicate marking, so the final read counts, EBV evidence, and
classification include all FASTQs for that sample. The final TSV includes an
`input_fastq_pairs` column recording how many pairs were combined.

For 10x single-cell Multiome ATAC data, supply the genomic `R1` and `R2`
FASTQs. This minimal screen does not use `I1`/`I2` and therefore reports
pseudobulk EBV evidence for the library, not EBV calls for individual cells.
For QNAMEs like the example ending in `:rCAGAC...`, leave
`barcode_aware_duplicates: true`: SAMtools then requires the eighth
colon-delimited field (the cell barcode) to match before marking two fragments
as duplicates. Set it to `false` for ordinary bulk ATAC read names without an
embedded cell barcode.

## Run

Create the minimal environment if needed:

```bash
mamba env create -f environment.yml
conda activate ebv-atac-minimal
```

Dry-run, then execute:

```bash
snakemake -n -p
snakemake --cores 8 --rerun-incomplete
```

If you do not activate the environment first, Snakemake can create the same
declared environment with `snakemake --use-conda --cores 8`.

For Slurm, request at least the CPU total configured for BWA plus sorting and
enough RAM for the hg38 BWA index and sort buffers (16 GB is a practical floor).

## Outputs

- `results/ebv_summary.tsv`: one row per sample with total read pairs, mapped
  reads, hg38 reads, EBV reads, conservative EBV fragments, RPM, control
  background, EBV breadth, and evidence classification. `ebv_covered_bases`
  is the union of bases covered by retained fragments, and
  `ebv_genome_breadth_fraction` divides that value by the EBV contig length.
  `ebv_rpm_over_negative_background` is sample EBV RPM divided by the maximum
  negative-control RPM. It is `Inf` for nonzero sample RPM over zero observed
  background and `NA` when no negative control was provided or when both values
  are zero.
- `results/evidence/SAMPLE.ebv_candidates.tsv`: every primary EBV candidate and
  its retain/reject reason, coordinates, MAPQ, CIGAR, aligned bases, edit
  distance, duplicate status, and optional repeat-region status.
- `results/bam/SAMPLE.hg38_EBV.bam[.bai]`: coordinate-sorted, duplicate-marked
  competitive alignment.
- `results/logs/`: BWA, sorting, duplicate, and validation logs.

## Conservative EBV definition

A retained EBV fragment must have exactly two primary records on the configured
EBV contig, be a proper pair, have both mates mapped with MAPQ at least 30,
have at least 30 aligned query bases per mate, pass the edit-distance limit,
and not be a marked duplicate. If `ebv_exclude_bed` is supplied, retained
fragments must also avoid those repeat or human-homologous intervals.

The workflow deliberately reports raw EBV mappings separately from retained
evidence. BWA may soft-clip adapters, but heavily clipped matches such as
`25S19M6S` fail the aligned-bases requirement.

`total_primary_reads` is counted from all primary BAM records. Because BWA
emits mapped and unmapped input records, it represents the FASTQ read count
after BWA has accepted the input; `total_read_pairs` is that value divided by
two. `human_mapped_reads` means mapped primary reads on any non-EBV contig, so
the combined reference should contain only hg38 (including its intended
decoys, if any) plus the named EBV contig.

Classification thresholds in `config/config.yaml` are screening defaults, not
universal biological truth. Matched negative controls establish background;
positive controls establish that the run had adequate sensitivity.
