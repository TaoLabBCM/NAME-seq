SAMPLES = ["PC9_W1_untreated"]
Read = ["R1", "R2"]
ref_genome = "/storage/taowu/home/xuwenl/ref_genome/human/gencode/GRCh38.primary_assembly.genome.fa"

rule all:
    input:
        expand("results_snp/trimmed/{sample}_{read}.fastq.gz", sample=SAMPLES, read=Read),
        expand("results_snp/bowtie/{sample}.sorted.bam", sample=SAMPLES),
        expand("results_snp/bcftools/{sample}.bed", sample=SAMPLES)



rule trimming:
    input:
        r1 = "raw_fastq/{sample}_R1.fastq.gz",
        r2 = "raw_fastq/{sample}_R2.fastq.gz"
    output:
        r1 = "results_snp/trimmed/{sample}_R1.fastq.gz",
        r2 = "results_snp/trimmed/{sample}_R2.fastq.gz"
        
    log:
        "logs_snp/cutadpt/{sample}.log"
    params:
        TRIM_OPTS = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -m 10"
    threads:
        32
    shell:
        """
        cutadapt {params.TRIM_OPTS} -j {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2}  >{log}
        """
        
rule bowtie2:
    input:
        "results_snp/trimmed/{sample}_R1.fastq.gz",
        "results_snp/trimmed/{sample}_R2.fastq.gz"
    params:
        index = "/storage/taowu/home/xuwenl/ref_genome/bowtie/indexes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
    output:
        temp("results_snp/bowtie/{sample}.sam")
    log:
        "logs_snp/bowtie/{sample}.log"
    threads: 30
    shell:
        "bowtie2 --local -x {params.index} -p {threads} -1 {input[0]} -2 {input[1]} -S {output} 2>{log}"

rule samtools:
    input:
        "results_snp/bowtie/{sample}.sam"
    output:
        "results_snp/bowtie/{sample}.sorted.bam"
    log:
        "logs_snp/samtools/{sample}.log"
    threads: 30
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2>{log}
        samtools index {output}
        """

rule bamCoverage:
    input: "results_snp/bowtie/{sample}.sorted.bam"
    output: "results_snp/deeptool/{sample}.bw"
    log: "logs_snp/deeptool/{sample}.log"
    shell:
        """
        bamCoverage -b {input} -o {output} -e -p 20 --ignoreDuplicates --binSize 30\
        --smoothLength 60 --effectiveGenomeSize 2494787188 \
        --normalizeUsing RPGC 2>{log}
        """
rule bcftools:
    input: "results_snp/bowtie/{sample}.sorted.bam"
    output: "results_snp/bcftools/{sample}.bed"
    log: "logs_snp/bcftools/{sample}.log"
    shell:
        """
        bcftools mpileup --threads 20 -f {ref_genome} {input} |\
        bcftools call --threads 20 -mv -Ob -o - |bcftools view --types snps -Ov | vcf2bed > {output}
        """

