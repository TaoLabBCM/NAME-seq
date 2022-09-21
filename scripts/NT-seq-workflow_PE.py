SAMPLES = ["HPJP26-Input-PCR", "HPJP26-Input-native"]
Read = ["R1", "R2"]
Strand = ["fwd", "rev"]
ref_genome = "../reference_genome/HPJP26.fasta"
ruleorder: fasta_AT_convert>trimming>dedupe>fastq_AT_convert>bowtie2_build>bowtie2_align>sam_to_original_pe>samtools_calmd>sam_tag_filtering>samtools_sort_index>igvtools
rule all:
    input:
        expand("results/AT_only_fasta/AT_only.fasta"),
        expand("bowtie2_index/index{ext}", ext=[".1.bt2", ".2.bt2", ".3.bt2",".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
        expand("results/trimmed_fastq/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=Read),
        expand("results/deduped_fastq/{sample}_{read}_dedupe.fastq", sample=SAMPLES, read=Read),
        expand("results/converted_fastq/{sample}_{read}_AT_only.fastq", sample=SAMPLES, read=Read),
        expand("results/converted_sam/{sample}_AT_only.sam", sample=SAMPLES),
        expand("results/original_sam/{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        expand("results/original_sam/{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        expand("results/original_sam/{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        expand("results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam",sample=SAMPLES, strand=Strand),
        expand("results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",sample=SAMPLES, strand=Strand),
        expand("results/readcount/{sample}_original_reads_{strand}_baq_filtered.wig",sample=SAMPLES, strand=Strand),
        expand("results/preprocessed_data/{sample}.csv",sample=SAMPLES),

rule fasta_AT_convert:
    input:
        ref_genome
    output:
        "results/AT_only_fasta/AT_only.fasta"
    shell:
        """
        python3 fasta_to_AT_only.py {input} {output}
        """

rule trimming:
    input:
        r1 = "../demo_data/{sample}_R1.fastq.gz",
        r2 = "../demo_data/{sample}_R2.fastq.gz",
    output:
        r1 = "results/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
    log:
        "logs/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 10"
    threads:
        32
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log}
        """
        
rule dedupe:
    input:
        r1 = "results/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed_fastq/{sample}_R2_trimmed.fastq.gz"
    output:
        r1 = "results/deduped_fastq/{sample}_R1_dedupe.fastq",
        r2 = "results/deduped_fastq/{sample}_R2_dedupe.fastq",
    log:
        "logs/dedupe/{sample}.log"
    shell:
        """
        clumpify.sh in1={input.r1} in2={input.r2} out={output.r1} out2={output.r2} dedupe 1>{log} 2>&1
        """

rule fastq_AT_convert:
    input:
        "results/deduped_fastq/{sample}_{read}_dedupe.fastq"
    output:
        "results/converted_fastq/{sample}_{read}_AT_only.fastq"
    shell:
        """
        python3 fastq_to_AT_only.py {input} {output}
        """

rule bowtie2_build:
    input:
        "results/AT_only_fasta/AT_only.fasta"
    output:
        multiext(
            "bowtie2_index/index",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log: "logs/bowtie2/index.log"
    shell:
        """
        bowtie2-build {input} bowtie2_index/index > {log} 2>&1
        """

rule bowtie2_align:
    input:
        r1 = "results/converted_fastq/{sample}_R1_AT_only.fastq",
        r2 = "results/converted_fastq/{sample}_R2_AT_only.fastq",
    output:
        "results/converted_sam/{sample}_AT_only.sam"
    log:
        "logs/bowtie/{sample}.log"
    threads: 30
    shell:
        "bowtie2 -x bowtie2_index/index -p {threads} -1 {input.r1} -2 {input.r2} -S {output} > {log} 2>&1"

rule sam_to_original_pe:
    input:
        sam = "results/converted_sam/{sample}_AT_only.sam",
        r1 = "results/deduped_fastq/{sample}_R1_dedupe.fastq",
        r2 = "results/deduped_fastq/{sample}_R2_dedupe.fastq",
    output:
        "results/original_sam/{sample}_original_reads_fwd.sam",
        "results/original_sam/{sample}_original_reads_rev.sam"
    shell:
        """
        python3 sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]} 
        """
rule samtools_calmd:
    input:
        "results/original_sam/{sample}_original_reads_{strand}.sam",
        ref_genome = ref_genome,
    output:
        "results/original_sam/{sample}_original_reads_{strand}.sorted.bam",
        "results/original_sam/{sample}_original_reads_{strand}_baq.sam"
    threads:
        10
    shell:
        """
        samtools sort {input[0]} > {output[0]}
        samtools calmd -@ {threads} -Ar {output[0]} {input.ref_genome} > {output[1]} 2>/dev/null
        """
rule sam_tag_filtering:
    input:
        fwd = "results/original_sam/{sample}_original_reads_fwd_baq.sam",
        rev = "results/original_sam/{sample}_original_reads_rev_baq.sam"
    output:
        fwd = "results/filtered_sam/{sample}_original_reads_fwd_baq_filtered.sam",
        rev = "results/filtered_sam/{sample}_original_reads_rev_baq_filtered.sam"

    shell:
        """
        python3 sam_tag_filtering.py {input.fwd} {output.fwd}
        python3 sam_tag_filtering.py {input.rev} {output.rev}
        """
rule samtools_sort_index:
    input:
        "results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam"
    output:
        "results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam"

    shell:
        """
        samtools sort -@ 20 {input} > {output[0]}
        samtools index {output[0]}
        """
rule igvtools:
    input:
        bam = "results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",
        ref_genome = ref_genome,
    output:
        "results/readcount/{sample}_original_reads_{strand}_baq_filtered.wig"
    shell:
        """
        igvtools count --bases -w 1 -e 0 {input.bam} {output} {input.ref_genome}
        """
rule data_proprocess:
    input:
        fwd = "results/readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = "results/readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        "results/preprocessed_data/{sample}.csv"
    shell:
        """
        python3 readcount_to_csv.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """
