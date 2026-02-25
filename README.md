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
