#!/usr/bin/env python3
#this line above tells the system to use Python 3 to execute this script
# when run from the terminal

import sys #the sys module allows us to access command-line arguments
# sys.argv is a list that stores arguments passed from the terminal
#sys.argv[0] = the script name
#sys.argv[1] = the first argument provided by the user (in this case, the FASTA file)

from Bio import SeqIO


THRESH = 1000
#only contigs longer than this threshold (1000 bp) will be counted

fasta = sys.argv[1]
#retrieves the FASTA filename provided from the terminal.
# Example usage:
# python script.py assembly.fasta
# In that case, sys.argv[1] = "assembly.fasta"

num_contigs = 0
#countr for contigs longer than THRESH

total_bp = 0
# Tracks total number of bps from contigs longer than THRESH

# go through the FASTA file record by record
for rec in SeqIO.parse(fasta, "fasta"):
    L = len(rec.seq)   #calculate length of current contig
    
    # only count contigs above the threshold length
    if L > THRESH:
        num_contigs += 1   
        total_bp += L  # add contig length to total base pairs


print(num_contigs)
print(total_bp)
