
# Overview

This repository contains a fully automated Snakemake pipeline for analyzing HCMV RNA-seq data from four SRA samples.

The pipeline performs:

1. Bowtie2 filtering (keeps mapped read pairs) to HCMV reference 
2. Read counts before and after filtering  
3. SPAdes assembly (k=99)  
4. Contig statistics (>1000 bp and total bp)  
5. Extraction of the longest contig  
6. BLAST of the longest contig against a locally built Betaherpesvirinae database  
7. Final summary written to:

Younis_PipelineReport.txt

Steps 2–5 are fully automated using Snakemake.

All file paths are relative so the pipeline can run on any system after cloning the repository.



# Repository Structure

Snakefile  
README.md  
scripts/  
    contig_stats.py  
    longest_contig.py  
data/  
    genome/  
        genome.fa  
    test_reads/  
        SRR5660030_1.fastq  
        SRR5660030_2.fastq  
        SRR5660033_1.fastq  
        SRR5660033_2.fastq  
        SRR5660044_1.fastq  
        SRR5660044_2.fastq  
        SRR5660045_1.fastq  
        SRR5660045_2.fastq  

The following folders are generated automatically when running the pipeline:

counts/
filtered_reads/
assemblies/
blast/
data/blastdb/
.snakemake/
*.bt2


# Required Software

You need the following installed:

- python3  
- snakemake  
- bowtie2  
- spades  
- blast+ (makeblastdb, blastn)  
- NCBI datasets CLI  
- unzip
- Biopython  



# Installing Software

Mac (Homebrew):

brew install snakemake  
brew install bowtie2  
brew install blast  
brew install spades  
brew install ncbi-datasets-cli  
brew install unzip  
(for biopython)
pip install biopython

OR Linux:

sudo apt install python3 python3-pip
pip3 install snakemake
sudo apt install bowtie2
sudo apt install spades
sudo apt install ncbi-blast+
sudo apt install unzip
pip3 install biopython




# Step 1: Download and Convert Reads (Manual Step)

The following SRA samples were retrieved:

SRR5660030 (Donor 1, 2dpi)  
SRR5660033 (Donor 1, 6dpi)  
SRR5660044 (Donor 3, 2dpi)  
SRR5660045 (Donor 3, 6dpi)  

Commands used:

mkdir -p data/full_reads data/test_reads  
cd data/full_reads  

fasterq-dump --split-files SRR5660030  
fasterq-dump --split-files SRR5660033  
fasterq-dump --split-files SRR5660044  
fasterq-dump --split-files SRR5660045  

cd ../..  



# Create Small Test FASTQs (Committed to GitHub)

To allow quicker testing small subsets were created below:

for s in SRR5660030 SRR5660033 SRR5660044 SRR5660045  
do  
  head -n 40000 data/full_reads/${s}_1.fastq > data/test_reads/${s}_1.fastq  
  head -n 40000 data/full_reads/${s}_2.fastq > data/test_reads/${s}_2.fastq  
done  

The pipeline defaults to using test_reads.



# Running the Pipeline (on the Test Data)
This single command executes the full workflow from Bowtie2 filtering through BLAST and making the final report.
From the main repo directory:

snakemake --cores 4 -p

This will:

- Build Bowtie2 index  
- Filter reads  
- Run SPAdes  
- Calculate contig stats  
- Download Betaherpesvirinae genomes automatically  
- Build local BLAST database  
- BLAST longest contig  
- Generate final report  

Output file:

Younis_PipelineReport.txt  

Runtime: < 2 minutes using test_reads.



# Running with Full Data

To run with full reads instead:

Open Snakefile and change:

READS_DIR = "data/test_reads"

to:

READS_DIR = "data/full_reads" (or whatever other filename desired)

Then run:

snakemake --cores 8 -p

(Could use 8 cores on a MacBook Air M4 for faster performance if you have a lot of RAM)



# Automated Betaherpesvirinae Database

The pipeline automatically downloads viral genomes using:

datasets download virus genome taxon 10357 --include genome

All FASTA files are concatenated into:

data/blastdb/betaherpes.fna

Then a local BLAST database is created using makeblastdb.




# Final Output

The final file Younis_PipelineReport.txt contains:

- Number of read pairs before and after Bowtie2 filtering  
- Number of contigs >1000 bp  
- Total bp in assembly  
- Top 5 BLAST hits for longest contig (best HSP only)  

This file is generated automatically by running Snakemake





