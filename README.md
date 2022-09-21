# NAME-seq
This pipeline is designed for data analysis of NAME-seq (Li, X., Guo, S., Cui, Y. et al. NAME-seq: a chemical-based sequencing method for genomic methylome profiling. Genome Biol 23, 122 (2022). https://doi.org/10.1186/s13059-022-02689-9).

## Requirements & Installation: 
Python3 \
Bowtie2 \
Cutadapt \
Snakemake \
Samtools \
Igvtools \
Biopython \
Docopt

Clone the repository: \
`git clone https://github.com/TaoLabBCM/NAME-seq`\
Add PATH of NAME-seq to your .bashrc file: \
`export PATH=$PATH:/PATH_OF_NAME-seq`\
Re log-in to terminal or use source command: \
`source ~/.bashrc`\
Create conda environment: \
`conda env create -f NAME-seq.yml`

## How to use it?
```
Usage:
    NTseq.py index [--index_name=index_name] [--num_threads=numthreads] (--input_fasta=fa_file)
    NTseq.py quant (--index=genome) (--input_fasta=fa_file) [--outpath=outpath] [--num_threads=numthreads] FILE...

Options:
    -h --help     Show this screen.
    --version     Show version.
```

## An example with demo data: 
Activate conda environment: \
`conda activate NAME-seq`\
Build AT only bowtie index:\
`NAME-seq.py index --index_name=HPJP26_AT_only --input_fasta=reference_genome/HPJP26.fasta`\
Run NAME-seq pipeline for single-end data:\
`NTseq.py quant --index=HPJP26_AT_only --input_fasta=reference_genome/HPJP26.fasta demo_data/HPJP26-NAME-seq_sampled_R1.fastq.gz`\
Run NAME-seq pipeline for paired-end data:\
`NTseq.py quant --index=HPJP26_AT_only --input_fasta=reference_genome/HPJP26.fasta demo_data/HPJP26-NAME-seq_sampled_R1.fastq.gz demo_data/HPJP26-NAME-seq_sampled_R2.fastq.gz`

### Expected output: 
CSV files containing base count at each position of reference genome for each sample, which can be used to calculate A to G ratio and C to T ratio for 6mA, 4mC, and 5mC motif calling.

### Test data: 
random sampled H. pylori JP26 NAME-seq data.

### Expected run time for demo
A few minutes using ten threads of Intel(R) Xeon(R) Gold 5120 CPU @ 2.20GHz

## Notes: 
To apply the script to other species, please build bowtie index with `NAME-seq.py index` command or copy existing bowtie index to index/ folder.
