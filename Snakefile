#https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-1-mapping-reads
#COMP383 Pipeline (Snakemake)
# - runs on TEST reads by default
# - generates Younis_PipelineReport.txt
#below is list of sample IDs we are processing.
#we use this so Snakemake can "expand" rules across all samples automatically
#instead of writing the same rule 4 times.
SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

#Where the input FASTQs live.
#switch this to data/full_reads for full data
READS_DIR = "data/test_reads"

#path to the HCMV reference genome fasta.
REF = "data/genome/genome.fa"
#taxonomy ID for Betaherpesvirinae on NCBI.
#we use a taxid because it’s one way to retrieve
# the right group automatically via the NCBI datasets tool.
BETAHERPES_TAXID = "10357"  # Betaherpesvirinae

# these are just file paths that we reuse in multiple rules.
# defining them once avoids hardcoding the same string everywhere.
BETAHERPES_ZIP = "data/blastdb/betaherpes_datasets.zip"
BETAHERPES_DIR = "data/blastdb/betaherpes_datasets"
BETAHERPES_FASTA = "data/blastdb/betaherpes.fna"

# BLAST DB prefix means BLAST will create multiple files like:
# betaherpes.nsq, betaherpes.nin, betaherpes.nhr, etc.
BLAST_DB_PREFIX = "data/blastdb/betaherpes"

# "rule all" defines the final targets we want Snakemake to produce
# When run snakemake it tries to build everything needed to create these.
rule all:
    input:
        "PipelineReport.txt"


# Step 2: Bowtie2 filtering + counts

rule bowtie2_index:
    input:
        REF
    output:
        REF + ".1.bt2"
    shell:
        #bowtie2-build output is actually multiple files, but one of them is
        # genome.fa.1.bt2. Snakemake just needs one output file to know the
        # index exists, so we use ".1.bt2" as the "flag" file.
        "bowtie2-build {input} {input}"

rule count_before:
    input:
        # NEW STYLE: lambda functions for inputs
        # Why this is used:
        # - Snakemake needs a way to turn the wc{sample} into real file paths
        # - Using lambda gives you full control over the filename pattern.
        # what it does:
        # - wc.sample is the wildcard value (example: "SRR5660030")
        # - returns "data/test_reads/SRR5660030_1.fastq" etc
        r1=lambda wc: f"{READS_DIR}/{wc.sample}_1.fastq",
        r2=lambda wc: f"{READS_DIR}/{wc.sample}_2.fastq"
    output:
        "counts/{sample}.before.txt"
    shell:
        # wc -l counts lines; FASTQ format is 4 lines per read
        # Since paired-end has read pairs, we count reads in R1 as the number of pairs
        "mkdir -p counts && "
        "echo $(( $(wc -l < {input.r1}) / 4 )) > {output}"

rule bowtie2_filter:
    input:
        r1=lambda wc: f"{READS_DIR}/{wc.sample}_1.fastq",
        r2=lambda wc: f"{READS_DIR}/{wc.sample}_2.fastq",
        idx=REF + ".1.bt2"
    output:
        r1_out="filtered_reads/{sample}_1.fastq",
        r2_out="filtered_reads/{sample}_2.fastq"
    shell:
        #triple-quoted multi-line shell string.
        # - makes long commands easier to read
        # - avoids unreadable one-line shells
        # bowtie2 --al-conc:
        # - writes ONLY the paired reads where BOTH mates align
        # - output uses a pattern with % meaning:
        #   filtered_reads/SRR5660030_1.fastq and ..._2.fastq
        # -S /dev/null means we throw away SAM output; we only want the filtered FASTQ pairs.
        """
        mkdir -p filtered_reads
        bowtie2 -x {REF} -1 {input.r1} -2 {input.r2} \
        --al-conc filtered_reads/{wildcards.sample}_%.fastq \
        -S /dev/null
        """

rule count_after:
    input:
        r1="filtered_reads/{sample}_1.fastq",
        r2="filtered_reads/{sample}_2.fastq"
    output:
        "counts/{sample}.after.txt"
    shell:
        #same logic as before: FASTQ = 4 lines per read.
        # Counting filtered R1 tells us how many *read pairs* remained.
        "echo $(( $(wc -l < {input.r1}) / 4 )) > {output}"


# Step 3: SPAdes assembly (k=99)

rule spades:
    input:
        r1="filtered_reads/{sample}_1.fastq",
        r2="filtered_reads/{sample}_2.fastq"
    output:
        "assemblies/{sample}/contigs.fasta"
    shell:
        # spades writes a whole directory of outputs
        # We point -o to the sample-specific folder
        # we redrect stdout+stderr into a log so the terminal isn't flooded.
        "mkdir -p assemblies/{wildcards.sample} && "
        "spades.py -1 {input.r1} -2 {input.r2} -k 99 -o assemblies/{wildcards.sample} "
        "> assemblies/{wildcards.sample}/spades.log 2>&1"


# Step 4: Contig stats (>1000 bp)

rule contig_stats:
    input:
        "assemblies/{sample}/contigs.fasta"
    output:
        "assemblies/{sample}/stats.txt"
    shell:
        #this runs the python script and saves output to a stats file
        "python3 scripts/contig_stats.py {input} > {output}"


# Step 5: Longest contig + BLAST
#Betaherpes FASTA using NCBI datasets

rule download_betaherpes_datasets:
    output:
        BETAHERPES_ZIP
    shell:
        #this uses the NCBI datasets CLI
        # "virus genome taxon 10357" = download viral genomes for Betaherpesvirinae.
        # --include genome = include the genome FASTA files
        # --filename = saves as a zip at the output path
        """
        mkdir -p data/blastdb
        datasets download virus genome taxon {BETAHERPES_TAXID} --include genome --filename {output}
        """

rule unzip_betaherpes_datasets:
    input:
        BETAHERPES_ZIP
    output:
        #directory()
        # - tells Snakemake "the output is a directory, not a single file"
        # - Snakemake then tracks the folder as the produced output
        directory(BETAHERPES_DIR)
    shell:
        #unzip the datasets zip into a clean folder.
        #rm -rf first prevents mixing old + new files
        """
        rm -rf {output}
        mkdir -p {output}
        unzip -q {input} -d {output}
        """

rule build_betaherpes_fasta:
    input:
        BETAHERPES_DIR
    output:
        BETAHERPES_FASTA
    shell:
        # We need ONE FASTA file to build a BLAST database.
        # datasets download produces many .fna files in subfolders.
        #NEW STYLE: find + xargs cat
        # - find ... -name "*.fna" lists all FASTA genome files
        # - sort makes the order consistent (reproducibility)
        # - xargs cat concatenates them all into one big .fna file
        """
        find {input} -type f -name "*.fna" -print0 | sort -z | xargs -0 cat > {output}
        """

rule make_blast_db:
    input:
        BETAHERPES_FASTA
    output:
        BLAST_DB_PREFIX + ".nsq"
    shell:
        # makeblastdb creates a nucleotide BLAST database from the combined FASTA
        # output is multiple files, but .nsq is one of them
        "makeblastdb -in {input} -out {BLAST_DB_PREFIX} -title betaherpes -dbtype nucl"

rule longest_contig:
    input:
        "assemblies/{sample}/contigs.fasta"
    output:
        "blast/{sample}.longest.fasta"
    shell:
        #Make blast folder then run script to extract the longest contig.
        "mkdir -p blast && python3 scripts/longest_contig.py {input} > {output}"

rule blast_top5:
    input:
        db=BLAST_DB_PREFIX + ".nsq",
        query="blast/{sample}.longest.fasta"
    output:
        "blast/{sample}.top5.tsv"
    shell:
        # blastn against our LOCAL Betaherpes DB
        # -max_target_seqs 5 = top 5 hits
        # -max_hsps 1 = keep best HSP only
        # -outfmt 6 = tabular output with the exact columns requested
        "blastn -query {input.query} -db {BLAST_DB_PREFIX} "
        "-max_hsps 1 -max_target_seqs 5 "
        "-outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' "
        "> {output}"


# Final: PipelineReport.txt

rule report:
    input:
        #expand()
        # Why:
        # - This automatically creates a list of filenames for each sample.
        # - Example: counts/SRR5660030.before.txt, counts/SRR5660033.before.txt, ...
        # Snakemake uses these lists to force all samples to finish before report runs
        before=expand("counts/{sample}.before.txt", sample=SAMPLES),
        after=expand("counts/{sample}.after.txt", sample=SAMPLES),
        stats=expand("assemblies/{sample}/stats.txt", sample=SAMPLES),
        top5=expand("blast/{sample}.top5.tsv", sample=SAMPLES)
    output:
        "PipelineReport.txt"
    shell:
        #builds the final report by looping over sample IDs in bash.
        # It writes:
        # 1) read pairs before/after
        # 2) contig count >1000 and total bp
        # 3) header + top 5 blast hits
        "echo '' > {output} && "
        "for s in " + " ".join(SAMPLES) + "; do "
        "b=$(cat counts/$s.before.txt); "
        "a=$(cat counts/$s.after.txt); "
        "echo \"Sample $s had $b read pairs before and $a read pairs after Bowtie2 filtering.\" >> {output}; "
        "num=$(head -n 1 assemblies/$s/stats.txt); "
        "bp=$(tail -n 1 assemblies/$s/stats.txt); "
        "echo \"In the assembly of sample $s, there are $num contigs > 1000 bp and $bp total bp.\" >> {output}; "
        "echo \"$s:\" >> {output}; "
        "echo -e \"sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\" >> {output}; "
        "cat blast/$s.top5.tsv >> {output}; "
        "echo '' >> {output}; "
        "done"
