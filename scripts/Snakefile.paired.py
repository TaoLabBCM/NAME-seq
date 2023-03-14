import sys
from os.path import join
INDEX = config["index"]
INPUT_DIR = config["input_path"]
OUTPUT_DIR = config["output_path"]
NUM_THREADS = int(config["num_threads"])
SAMPLES = config["basename_list"]
ref_genome = config["fasta"]

Read = ["R1", "R2"]
Strand = ["fwd", "rev"]

# ruleorder: dedupe>trimming>fastq_AT_convert>bowtie2_align>mapped_fastq_filtering>sam_to_original_pe>samtools_calmd>sam_tag_filtering>samtools_sort_index>igvtools
rule all:
    input:
        # expand(OUTPUT_DIR + "deduped_fastq/{sample}_{read}_dedupe.fastq.gz", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed2.fastq", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "converted_fastq/{sample}_{read}_trimmed2_AT_only.fastq", sample=SAMPLES, read=Read),
        # expand(OUTPUT_DIR + "converted_sam/{sample}_trimmed2_AT_only.sam", sample=SAMPLES),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam",sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",sample=SAMPLES, strand=Strand),
        # expand(OUTPUT_DIR + "readcount/{sample}_original_reads_{strand}_baq_filtered.wig",sample=SAMPLES, strand=Strand),
        expand(OUTPUT_DIR + "preprocessed_data/{sample}_adenine.csv",sample=SAMPLES),
        expand(OUTPUT_DIR + "preprocessed_data/{sample}_cytosine.csv",sample=SAMPLES),

rule dedupe:
    input:
        r1 = "raw_fastq/{sample}_R1.fastq.gz",
        r2 = "raw_fastq/{sample}_R2.fastq.gz",
    output:
        r1 = temp(OUTPUT_DIR + "deduped_fastq/{sample}_R1_dedupe.fastq.gz"),
        r2 = temp(OUTPUT_DIR + "deduped_fastq/{sample}_R2_dedupe.fastq.gz"),
    log:
        "logs/dedupe/{sample}.log"
    shell:
        """
        clumpify.sh -Xmx20G in1={input.r1} in2={input.r2} out={output.r1} out2={output.r2} dedupe 1>{log} 2>&1
        """
# rule pre_trimming:
#     input:
#         r1 = OUTPUT_DIR + "deduped_fastq/{sample}_R1_dedupe.fastq.gz",
#         r2 = OUTPUT_DIR + "deduped_fastq/{sample}_R2_dedupe.fastq.gz",
#     output:
#         r1 = OUTPUT_DIR + "trimmed_fastq/{sample}_R1_pretrimmed.fastq.gz",
#         r2 = OUTPUT_DIR + "trimmed_fastq/{sample}_R2_pretrimmed.fastq.gz",
#     log:
#         "logs/cutadpt/{sample}.log"
#     params:
#         TRIM_OPTS = "-m 10 -u -75 -U -75"
#     threads:
#         32
#     shell:
#         """
#         cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log}
#         """
        
rule trimming:
    input:
        r1 = OUTPUT_DIR + "deduped_fastq/{sample}_R1_dedupe.fastq.gz",
        r2 = OUTPUT_DIR + "deduped_fastq/{sample}_R2_dedupe.fastq.gz",
    output:
        r1 = temp(OUTPUT_DIR + "trimmed_fastq/{sample}_R1_trimmed.fastq.gz"),
        r2 = temp(OUTPUT_DIR + "trimmed_fastq/{sample}_R2_trimmed.fastq.gz"),
    log:
        "logs/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 10"
    threads: NUM_THREADS
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} -u 9 -U 9 1>{log}
        """

rule trimming2:
    input:
        r1 = OUTPUT_DIR + "trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        r2 = OUTPUT_DIR + "trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
    output:
        r1 = temp(OUTPUT_DIR + "trimmed_fastq/{sample}_R1_trimmed2.fastq"),
        r2 = temp(OUTPUT_DIR + "trimmed_fastq/{sample}_R2_trimmed2.fastq"),
    log:
        "logs/cutadpt2/{sample}.log"
    params:
        TRIM_OPTS = "-m 10 -u -9 -U -9"
    threads: NUM_THREADS
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log}
        """

rule fastq_AT_convert:
    input:
        OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed2.fastq"
    output:
        temp(OUTPUT_DIR + "converted_fastq/{sample}_{read}_trimmed2_AT_only.fastq")
    shell:
        """
        python3 fastq_to_AT_only.py {input} {output}
        """

rule bowtie2_align:
    input:
        r1 = OUTPUT_DIR + "converted_fastq/{sample}_R1_trimmed2_AT_only.fastq",
        r2 = OUTPUT_DIR + "converted_fastq/{sample}_R2_trimmed2_AT_only.fastq",
    output:
        temp(OUTPUT_DIR + "converted_sam/{sample}_trimmed2_AT_only.sam")
    params:
        index = INDEX,
    log:
        "logs/bowtie/{sample}.log"
    threads: NUM_THREADS
    shell:
        "bowtie2 -x {params} -p {threads} -1 {input.r1} -2 {input.r2} -S {output} > {log} 2>&1"

rule mapped_fastq_filtering:
    input:
        sam = OUTPUT_DIR + "converted_sam/{sample}_trimmed2_AT_only.sam",
        fastq = OUTPUT_DIR + "trimmed_fastq/{sample}_{read}_trimmed2.fastq"
    output:
        fastq = temp(OUTPUT_DIR + "mapped_fastq/{sample}_{read}_mapped.fastq"),
        readid = temp('results/mapped_fastq/{sample}_{read}_mapped.readid')
    shell:
        """
        samtools view -F 4 {input.sam} | cut -f 1 | sort | uniq > {output.readid}
        filterbyname.sh -Xmx20G in={input.fastq} out={output.fastq} names={output.readid} include=t
        """

rule sam_to_original_pe:
    input:
        sam = OUTPUT_DIR + "converted_sam/{sample}_trimmed2_AT_only.sam",
        r1 = OUTPUT_DIR + "mapped_fastq/{sample}_R1_mapped.fastq",
        r2 = OUTPUT_DIR + "mapped_fastq/{sample}_R2_mapped.fastq",
    output:
        temp(OUTPUT_DIR + "original_sam/{sample}_original_reads_fwd.sam"),
        temp(OUTPUT_DIR + "original_sam/{sample}_original_reads_rev.sam")
    shell:
        """
        python3 sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]} 
        """
rule samtools_calmd:
    input:
        OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sam",
        ref_genome = ref_genome,
    output:
        temp(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}.sorted.bam"),
        temp(OUTPUT_DIR + "original_sam/{sample}_original_reads_{strand}_baq.sam")
    threads:
        10
    shell:
        """
        samtools sort {input[0]} > {output[0]}
        samtools calmd -@ {threads} {output[0]} {input.ref_genome} > {output[1]} 2>/dev/null
        """
rule sam_tag_filtering:
    input:
        fwd = OUTPUT_DIR + "original_sam/{sample}_original_reads_fwd_baq.sam",
        rev = OUTPUT_DIR + "original_sam/{sample}_original_reads_rev_baq.sam"
    output:
        fwd = temp(OUTPUT_DIR + "filtered_sam/{sample}_original_reads_fwd_baq_filtered.sam"),
        rev = temp(OUTPUT_DIR + "filtered_sam/{sample}_original_reads_rev_baq_filtered.sam")

    shell:
        """
        python3 sam_tag_filtering.py {input.fwd} {output.fwd}
        python3 sam_tag_filtering.py {input.rev} {output.rev}
        """
rule samtools_sort_index:
    input:
        OUTPUT_DIR + "filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam"
    output:
        temp(OUTPUT_DIR + "filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam")

    shell:
        """
        samtools sort -@ 20 {input} > {output[0]}
        samtools index {output[0]}
        """
rule igvtools:
    input:
        bam = OUTPUT_DIR + "filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",
        ref_genome = ref_genome,
    output:
        temp(OUTPUT_DIR + "readcount/{sample}_original_reads_{strand}_baq_filtered.wig")
    shell:
        """
        igvtools count --bases -w 1 -e 0 {input.bam} {output} {input.ref_genome} --includeDuplicates
        """
rule data_proprocess:
    input:
        fwd = OUTPUT_DIR + "readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = OUTPUT_DIR + "readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        OUTPUT_DIR + "preprocessed_data/{sample}_adenine.csv"
    shell:
        """
        python3 readcount_to_csv.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """
        
rule data_proprocess_cytosine:
    input:
        fwd = OUTPUT_DIR + "readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = OUTPUT_DIR + "readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        OUTPUT_DIR + "preprocessed_data/{sample}_cytosine.csv"
    shell:
        """
        python3 readcount_to_csv-cytosine.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """