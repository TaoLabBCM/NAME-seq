SAMPLES = [
           "PC9_W1_NAME-seq",
          ]
Read = ["R1", "R2"]
Strand = ["fwd", "rev"]
ref_genome = "/storage/taowu/home/xuwenl/ref_genome/human/gencode/GRCh38.primary_assembly.genome.fa"
bowtie_index = "/storage/taowu/home/xuwenl/ref_genome/bowtie/indexes/GRCh38_AT_only/GRCh38_AT_only"
# ruleorder: dedupe>trimming>trimming2>fastq_AT_convert>bowtie2_align>mapped_fastq_filtering>sam_to_original_pe>samtools_calmd>sam_tag_filtering>samtools_sort_index>igvtools
rule all:
    input:
        expand("results/deduped_fastq/{sample}_{read}_dedupe.fastq.gz", sample=SAMPLES, read=Read),
        expand("results/trimmed_fastq/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=Read),
        # expand("results/trimmed_fastq/{sample}_{read}_trimmed2.fastq", sample=SAMPLES, read=Read),
        # expand("results/converted_fastq/{sample}_{read}_trimmed2_AT_only.fastq", sample=SAMPLES, read=Read),
        # expand("results/converted_sam/{sample}_trimmed2_AT_only.sam", sample=SAMPLES),
        # expand("results/original_sam/{sample}_original_reads_{strand}.sam", sample=SAMPLES, strand=Strand),
        # expand("results/original_sam/{sample}_original_reads_{strand}.sorted.bam", sample=SAMPLES, strand=Strand),
        # expand("results/original_sam/{sample}_original_reads_{strand}_baq.sam", sample=SAMPLES, strand=Strand),
        # expand("results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sam",sample=SAMPLES, strand=Strand),
        expand("results/filtered_sam/{sample}_original_reads_{strand}_baq_filtered.sorted.bam",sample=SAMPLES, strand=Strand),
        expand("results/readcount/{sample}_original_reads_{strand}_baq_filtered.wig",sample=SAMPLES, strand=Strand),
        expand("results/preprocessed_data/{sample}_adenine.csv",sample=SAMPLES),
        expand("results/preprocessed_data/{sample}_cytosine.csv",sample=SAMPLES),

rule dedupe:
    input:
        r1 = "raw_fastq/{sample}_R1.fastq.gz",
        r2 = "raw_fastq/{sample}_R2.fastq.gz",
    output:
        r1 = temp("results/deduped_fastq/{sample}_R1_dedupe.fastq.gz"),
        r2 = temp("results/deduped_fastq/{sample}_R2_dedupe.fastq.gz"),
    log:
        "logs/dedupe/{sample}.log"
    shell:
        """
        clumpify.sh -Xmx20G in1={input.r1} in2={input.r2} out={output.r1} out2={output.r2} dedupe 1>{log} 2>&1
        """
# rule pre_trimming:
#     input:
#         r1 = "results/deduped_fastq/{sample}_R1_dedupe.fastq.gz",
#         r2 = "results/deduped_fastq/{sample}_R2_dedupe.fastq.gz",
#     output:
#         r1 = "results/trimmed_fastq/{sample}_R1_pretrimmed.fastq.gz",
#         r2 = "results/trimmed_fastq/{sample}_R2_pretrimmed.fastq.gz",
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
        r1 = "results/deduped_fastq/{sample}_R1_dedupe.fastq.gz",
        r2 = "results/deduped_fastq/{sample}_R2_dedupe.fastq.gz",
    output:
        r1 = temp("results/trimmed_fastq/{sample}_R1_trimmed.fastq.gz"),
        r2 = temp("results/trimmed_fastq/{sample}_R2_trimmed.fastq.gz"),
    log:
        "logs/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 10"
    threads:
        32
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} -u 9 -U 9 1>{log}
        """

rule trimming2:
    input:
        r1 = "results/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
    output:
        r1 = temp("results/trimmed_fastq/{sample}_R1_trimmed2.fastq"),
        r2 = temp("results/trimmed_fastq/{sample}_R2_trimmed2.fastq"),
    log:
        "logs/cutadpt2/{sample}.log"
    params:
        TRIM_OPTS = "-m 10 -u -66 -U -66"
    threads:
        32
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} 1>{log}
        """

rule fastq_AT_convert:
    input:
        "results/trimmed_fastq/{sample}_{read}_trimmed2.fastq"
    output:
        temp("results/converted_fastq/{sample}_{read}_trimmed2_AT_only.fastq")
    shell:
        """
        python3 fastq_to_AT_only.py {input} {output}
        """

rule bowtie2_align:
    input:
        r1 = "results/converted_fastq/{sample}_R1_trimmed2_AT_only.fastq",
        r2 = "results/converted_fastq/{sample}_R2_trimmed2_AT_only.fastq",
    output:
        temp("results/converted_sam/{sample}_trimmed2_AT_only.sam")
    params:
        bowtie_index
    log:
        "logs/bowtie/{sample}.log"
    threads: 30
    shell:
        "bowtie2 -x {params} -p {threads} -1 {input.r1} -2 {input.r2} -S {output} > {log} 2>&1"

rule mapped_fastq_filtering:
    input:
        sam = "results/converted_sam/{sample}_trimmed2_AT_only.sam",
        fastq = "results/trimmed_fastq/{sample}_{read}_trimmed2.fastq"
    output:
        fastq = temp("results/mapped_fastq/{sample}_{read}_mapped.fastq"),
        readid = temp('results/mapped_fastq/{sample}_{read}_mapped.readid')
    shell:
        """
        samtools view -F 4 {input.sam} | cut -f 1 | sort | uniq > {output.readid}
        filterbyname.sh -Xmx20G in={input.fastq} out={output.fastq} names={output.readid} include=t
        """

rule sam_to_original_pe:
    input:
        sam = "results/converted_sam/{sample}_trimmed2_AT_only.sam",
        r1 = "results/mapped_fastq/{sample}_R1_mapped.fastq",
        r2 = "results/mapped_fastq/{sample}_R2_mapped.fastq",
    output:
        temp("results/original_sam/{sample}_original_reads_fwd.sam"),
        temp("results/original_sam/{sample}_original_reads_rev.sam")
    shell:
        """
        python3 sam_to_original_pe.py {input.sam} {input.r1} {input.r2} {output[0]} {output[1]} 
        """
rule samtools_calmd:
    input:
        "results/original_sam/{sample}_original_reads_{strand}.sam",
        ref_genome = ref_genome,
    output:
        temp("results/original_sam/{sample}_original_reads_{strand}.sorted.bam"),
        temp("results/original_sam/{sample}_original_reads_{strand}_baq.sam")
    threads:
        10
    shell:
        """
        samtools sort {input[0]} > {output[0]}
        samtools calmd -@ {threads} {output[0]} {input.ref_genome} > {output[1]} 2>/dev/null
        """
rule sam_tag_filtering:
    input:
        fwd = "results/original_sam/{sample}_original_reads_fwd_baq.sam",
        rev = "results/original_sam/{sample}_original_reads_rev_baq.sam"
    output:
        fwd = temp("results/filtered_sam/{sample}_original_reads_fwd_baq_filtered.sam"),
        rev = temp("results/filtered_sam/{sample}_original_reads_rev_baq_filtered.sam")

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
        igvtools count --bases -w 1 -e 0 {input.bam} {output} {input.ref_genome} --includeDuplicates
        """
rule data_proprocess:
    input:
        fwd = "results/readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = "results/readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        "results/preprocessed_data/{sample}_adenine.csv"
    shell:
        """
        python3 readcount_to_csv.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """
        
rule data_proprocess_cytosine:
    input:
        fwd = "results/readcount/{sample}_original_reads_fwd_baq_filtered.wig",
        rev = "results/readcount/{sample}_original_reads_rev_baq_filtered.wig",
        ref_genome = ref_genome,
    output:
        "results/preprocessed_data/{sample}_cytosine.csv"
    shell:
        """
        python3 readcount_to_csv-cytosine.py {input.fwd} {input.rev} {input.ref_genome} {output}
        """