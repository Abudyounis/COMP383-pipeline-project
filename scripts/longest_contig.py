#!/usr/bin/env python3
# uses the system's Python 3 when run from terminal (automatically generated in Spyder)

import sys
#the sys module allows access to:
#command-line arguments (sys.argv)
#standard output (sys.stdout)
#exit functions (sys.exit)
from Bio import SeqIO

fasta = sys.argv[1]
#retrieves the FASTA filename passed from the terminal.
# example:
# python longest_contig.py assembly.fasta
# sys.argv[1] = "assembly.fasta"

best = None
# gunna store the longest contig record found

best_len = -1
#tracks length of longest contig seen so far
# initialized to -1 so that any real contig length (>= 0) will replace it

#go through each sequence record in the FASTA file
for rec in SeqIO.parse(fasta, "fasta"):
    L = len(rec.seq)  #calculate length of current contig
    
    # if current contig is longer than the previous best
    if L > best_len:
        best_len = L
        best = rec   # store this contig as the new longest

# if no contigs were found in the file
if best is None:
    sys.exit(1)
    # exit with status code 1 (shows error to the shell / Snakemake)

#write the longest contig to standard output in FASTA format
SeqIO.write(best, sys.stdout, "fasta")
