import sys
from os.path import join
INDEX = config["index"]
INPUT_DIR = config["input_path"]
OUTPUT_DIR = config["output_path"]
NUM_THREADS = config["num_threads"]
SAMPLES = config["basename_list"]
ref_genome = config["fasta"]

Strand = ["fwd", "rev"]

ruleorder: dedupe>trimming>fastq_AT_convert>bowtie2_align>sam_to_original_se>samtools_calmd>sam_tag_filtering>samtools_sort_index>igvtools
rule all:
    input:
        expand("NTseq_output/trimmed_fastq/{sample}_trimmed.fastq", sample=SAMPLES),
        expand("NTseq_output/deduped_fastq/{sample}_dedupe.fastq.gz", sample=SAMPLES),
        expand("NTseq_output/converted_fastq/{sample}_AT_only.fastq", sample=SAMPLES),
        expand("NTseq_output/converted_sam/{sample}_AT_only.sam", sample=SAMPLES),
        expand("NTseq_output/original_sam/{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        expand("NTseq_output/original_sam/{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        expand("NTseq_output/original_sam/{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        expand("NTseq_output/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam",sample=SAMPLES, strand=Strand),
        expand("NTseq_output/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",sample=SAMPLES, strand=Strand),
        expand("NTseq_output/readcount/{sample}_original_reads_{strand}_baq_filtered.wig",sample=SAMPLES, strand=Strand),
        expand("NTseq_output/preprocessed_data/{sample}.csv",sample=SAMPLES)


rule dedupe:
    input:
        r1 = "demo_data/{sample}.fastq.gz"
    output:
        r1 = "NTseq_output/deduped_fastq/{sample}_dedupe.fastq.gz"
    log:
        "NTseq_output/logs/dedupe/{sample}.log"
    shell:
        """
        clumpify.sh in={input.r1} out={output.r1} dedupe 1>{log} 2>&1
        """

rule trimming:
    input:
        r1 = "NTseq_output/deduped_fastq/{sample}_dedupe.fastq.gz"
    output:
        r1 = "NTseq_output/trimmed_fastq/{sample}_trimmed.fastq"
    log:
        "NTseq_output/logs/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 10"
    threads:
        NUM_THREADS
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} {input.r1} 1>{log}
        """
        

rule fastq_AT_convert:
    input:
        "NTseq_output/trimmed_fastq/{sample}_trimmed.fastq"
    output:
        "NTseq_output/converted_fastq/{sample}_AT_only.fastq"
    shell:
        """
        python scripts/fastq_to_AT_only.py {input} {output}
        """

rule bowtie2_align:
    input:
        r1 = "NTseq_output/converted_fastq/{sample}_AT_only.fastq"
    params:
        index = INDEX,
    output:
        "NTseq_output/converted_sam/{sample}_AT_only.sam"
    log:
        "NTseq_output/logs/bowtie/{sample}.log"
    threads: NUM_THREADS
    shell:
        "bowtie2 -x {params.index} -p {threads} -U {input.r1} -S {output} > {log} 2>&1"

rule sam_to_original_se:
    input:
        sam = "NTseq_output/converted_sam/{sample}_AT_only.sam",
        r1 = "NTseq_output/trimmed_fastq/{sample}_trimmed.fastq"
    output:
        "NTseq_output/original_sam/{sample}_original_reads_fwd.sam",
        "NTseq_output/original_sam/{sample}_original_reads_rev.sam"
    shell:
        """
        python scripts/sam_to_original_se.py {input.sam} {input.r1} {output[0]} {output[1]} 
        """
rule samtools_calmd:
    input:
        "NTseq_output/original_sam/{sample}_original_reads_{strand}.sam",
        ref_genome = ref_genome,
    output:
        "NTseq_output/original_sam/{sample}_original_reads_{strand}.sorted.bam",
        "NTseq_output/original_sam/{sample}_original_reads_{strand}_baq.sam"
    threads:
        NUM_THREADS
    shell:
        """
        samtools sort {input[0]} > {output[0]}
        samtools calmd -@ {threads} -Ar {output[0]} {input.ref_genome} > {output[1]} 2>/dev/null
        """
rule sam_tag_filtering:
    input:
        fwd = "NTseq_output/original_sam/{sample}_original_reads_fwd_baq.sam",
        rev = "NTseq_output/original_sam/{sample}_original_reads_rev_baq.sam"
    output:
        fwd = "NTseq_output/filtered_sam/{sample}_original_reads_fwd_baq_filtered.sam",
        rev = "NTseq_output/filtered_sam/{sample}_original_reads_rev_baq_filtered.sam"

    shell:
        """
        python scripts/sam_tag_filtering.py {input.fwd} {output.fwd}
        python scripts/sam_tag_filtering.py {input.rev} {output.rev}
        """
rule samtools_sort_index:
    input:
        "NTseq_output/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam"
    output:
        "NTseq_output/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam"

    shell:
        """
        samtools sort -@ 20 {input} > {output[0]}
        samtools index {output[0]}
        """
rule igvtools:
    input:
        bam = "NTseq_output/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",
        ref_genome = ref_genome,
    output:
        "NTseq_output/readcount/{sample}_original_reads_{strand}_baq_filtered.wig"
    shell:
        """
        igvtools count --bases -w 1 -e 0 {input.bam} {output} {input.ref_genome}
        """
rule data_proprocess:
    input:
        fwd = "NTseq_output/readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = "NTseq_output/readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        "NTseq_output/preprocessed_data/{sample}.csv"
    shell:
        """
        python scripts/readcount_to_csv.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """
