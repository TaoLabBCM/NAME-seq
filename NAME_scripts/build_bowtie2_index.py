ref_genome = "../reference_genome/HPJP26.fasta"
rule all:
    input:
        expand("results/AT_only_fasta/AT_only.fasta"),
        expand("bowtie2_index/index{ext}", ext=[".1.bt2", ".2.bt2", ".3.bt2",".4.bt2", ".rev.1.bt2", ".rev.2.bt2"])

rule fasta_AT_convert:
    input:
        ref_genome
    output:
        "results/AT_only_fasta/AT_only.fasta"
    shell:
        """
        python3 fasta_to_AT_only.py {input} {output}
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
