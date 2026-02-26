SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
READS_DIR = "data/test_reads"
REF = "data/genome/genome.fa"

rule all:
    input:
        "Younis_PipelineReport.txt"

# Step 2: Bowtie2 filtering + counts


rule bowtie2_index:
    input:
        REF
    output:
        REF + ".1.bt2"
    shell:
        "bowtie2-build {input} {input}"

rule count_before:
    input:
        r1=lambda wc: f"{READS_DIR}/{wc.sample}_1.fastq",
        r2=lambda wc: f"{READS_DIR}/{wc.sample}_2.fastq"
    output:
        "counts/{sample}.before.txt"
    shell:
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
        "echo $(( $(wc -l < {input.r1}) / 4 )) > {output}"


# Step 3: SPAdes assembly (k=99)


rule spades:
    input:
        r1="filtered_reads/{sample}_1.fastq",
        r2="filtered_reads/{sample}_2.fastq"
    output:
        "assemblies/{sample}/contigs.fasta"
    shell:
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
        "python3 scripts/contig_stats.py {input} > {output}"

# Step 5: Longest contig + BLAST (top 5, best HSP only)

rule make_blast_db:
    input:
        "data/blastdb/betaherpes.fna"
    output:
        "data/blastdb/betaherpes.nsq"
    shell:
        "makeblastdb -in {input} -out data/blastdb/betaherpes -title betaherpes -dbtype nucl"

rule longest_contig:
    input:
        "assemblies/{sample}/contigs.fasta"
    output:
        "blast/{sample}.longest.fasta"
    shell:
        "mkdir -p blast && python3 scripts/longest_contig.py {input} > {output}"

rule blast_top5:
    input:
        db="data/blastdb/betaherpes.nsq",
        query="blast/{sample}.longest.fasta"
    output:
        "blast/{sample}.top5.tsv"
    shell:
        "blastn -query {input.query} -db data/blastdb/betaherpes "
        "-max_hsps 1 -max_target_seqs 5 "
        "-outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' "
        "> {output}"

# Final: PipelineReport.txt


rule report:
    input:
        before=expand("counts/{sample}.before.txt", sample=SAMPLES),
        after=expand("counts/{sample}.after.txt", sample=SAMPLES),
        stats=expand("assemblies/{sample}/stats.txt", sample=SAMPLES),
        top5=expand("blast/{sample}.top5.tsv", sample=SAMPLES)
    output:
        "Younis_PipelineReport.txt"
    shell:
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