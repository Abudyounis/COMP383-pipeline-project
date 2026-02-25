# COMP383 Python Snakemake Pipeline Project

# Overview
This repo contains a Snakemake pipeline that automates COMP383 Pipeline Project steps:
- Bowtie2 mapping/filtering to HCMV reference
- SPAdes assembly (k=99)
- Assembly contig stats (>1000 bp) and total bp
- Longest contig BLAST against a local Betaherpesvirinae database
- Writes results to `results/reports/PipelineReport.txt`

# Required tools:
- snakemake
- python3
- bowtie2
- spades.py
- blast+ (makeblastdb, blastn)
- biopython
# Project structure
This repo follows the Snakemake tutorial directory layout:

- `Snakefile` — Snakemake workflow (rules)
- `data/` — input data
  - `data/genome.fa` — reference genome FASTA used by bwa
  - `data/samples/` — FASTQ inputs named `{sample}.fastq`
- `mapped_reads/` — BAM outputs from bwa mapping (`{sample}.bam`)
- `sorted_reads/` — sorted BAM outputs and BAM indexes (`{sample}.bam`, `{sample}.bam.bai`)
- `scripts/` — Python scripts used by Snakemake rules

# Step 1: Download and convert reads (manual)
We downloaded and converted the following SRA runs to paired-end FASTQ using `fasterq-dump --split-files`:

- SRR5660030 (Donor 1, 2dpi)
- SRR5660033 (Donor 1, 6dpi)
- SRR5660044 (Donor 3, 2dpi)
- SRR5660045 (Donor 3, 6dpi)

Commands used:
mkdir -p data/full_reads data/test_reads
cd data/full_reads
fasterq-dump --split-files SRR5660030
fasterq-dump --split-files SRR5660033
fasterq-dump --split-files SRR5660044
fasterq-dump --split-files SRR5660045
cd ../..

# create small test FASTQs (committed to GitHub)
for s in SRR5660030 SRR5660033 SRR5660044 SRR5660045
do
  head -n 40000 data/full_reads/${s}_1.fastq > data/test_reads/${s}_1.fastq
  head -n 40000 data/full_reads/${s}_2.fastq > data/test_reads/${s}_2.fastq
done
