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
