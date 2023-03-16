#!/usr/bin/env python
"""NAME-seq: a sequencing method for quantitative mapping of DNA amino methylation (N6-methyladenine and N4-methylcytosine) at single-base resolution

Usage:
    NAMEseq.py index [--index_name=index_name] [--num_threads=numthreads] (--input_fasta=fa_file)
    NAMEseq.py quant (--index=genome) (--input_fasta=fa_file) [--outpath=outpath] [--num_threads=numthreads] FILE...

Options:
    -h --help     Show this screen.
    --version     Show version.
    --index_name  Name for the Bowtie2 index output
    -t, --num_threads Num of threads to be used
    --input_fasta Input reference genome file
    --index       Built Bowtie2 index location
    -o, --outpath     Path to store the output files

"""

from docopt import docopt
from glob import glob
from itertools import combinations
import logging
import sys
import os
import tempfile
import gzip

def is_fastq(file_name):
    file_name = file_name.lower()
    if file_name.endswith("fastq"): return True
    if file_name.endswith("fq"): return True
    if file_name.endswith("fq.gz"): return True
    if file_name.endswith("fastq.gz"): return True
    return False

def get_basename_noext(file_name):
    # return os.path.basename(file_name.split('.')[0])
    return os.path.basename(file_name).split('.')[0]


def get_ext(file_name):
    # return ".".join(os.path.basename(file_name).split('.')[1:])
    return ".".join(os.path.basename(file_name).split('.')[1:])

def longest_prefix(a, b):
    a, b = (b,a) if len(a) > len(b) else (a,b)
    return b[:max([i for i in range(len(a),-1,-1) if b.startswith(a[:i])]+[0])]

def get_first_readid(file_name):
    if file_name.lower().endswith(".gz"):
        with gzip.open(file_name, "rb") as inp:
            return inp.readline().decode("utf8").replace("/"," ").split()[0]
    else:
        with open(file_name, "r") as inp:
            return inp.readline().replace("/"," ").split()[0]


def correct_ext(ext):
    return "fastq.gz" if "gz" in ext else "fastq"

def collect_FASTQ_files(FILE):
    print(FILE)
    fastq_files = set()
    for file in FILE:
        collected_files = glob(file)
        if len(collected_files) == 0:
            logging.error("Failed to read : '{}' is not existed.".format(file))
            sys.exit(1)

        elif len(collected_files) == 1 and not is_fastq(collected_files[0]):
            logging.info("This script assumes that '{}' is a directory, and will search any FASTQ file in the directory.".format(collected_files[0]))
            collected_files = glob(collected_files[0]+"/*")

        black_list = set()
        for col_file in collected_files:
            if not is_fastq(col_file):
                # logging.error("Failed to read : '{}' seems not to be a FASTQ file.".format(file))
                logging.warn("NAMEseq found '{}' file, but this seems not to be a FASTQ file.".format(col_file))
                black_list.add(col_file)

        fastq_files |= set([os.path.abspath(file) for file in collected_files if file not in black_list])

    tmp_dir = tempfile.mkdtemp()

    paired = dict([ (file,set()) for file in fastq_files ])
    for file_a, file_b in combinations(fastq_files, 2):
        """
        if get_ext(file_a) != get_ext(file_b):
            logging.info("File extensions of all files must be same.")
            sys.exit(1)
        """

        if get_first_readid(file_a) == get_first_readid(file_b):
            paired[file_a].add(file_b)
            paired[file_b].add(file_a)

    is_paired = True
    for file in paired:
        if len(paired[file]) == 0: is_paired = False



    for file in paired:
        if is_paired and len(paired[file]) > 1:
            logging.error("One of input files can be paired to multiple files")
            sys.exit(1)
        elif not is_paired and len(paired[file]) > 0:
            logging.error("A paired-end sample and a single-end sample are placed together.")
            sys.exit(1)

    file_list = []
    basename_list = []
    if is_paired:
        for file, matched in paired.items():
            if len(matched) == 0:
                logging.error("Input dataset is supposed as a paired-ends dataset, but some files do not have thier pair.")
                sys.exit(1)

        logging.info("The input dataset is considered as a paired-ends dataset.")
        for file in sorted(paired):
            a = file
            b = list(paired[file])[0]
            print(paired)
            print("a:{}".format(a))
            print("b:{}".format(b))
            if a > b: continue
            trim_a = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R1.{}".format(correct_ext(os.path.basename(a).split('.')[1:]))
            trim_b = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R2.{}".format(correct_ext(os.path.basename(b).split('.')[1:]))
            print(tmp_dir)
            os.symlink(os.path.abspath(a), os.path.join(tmp_dir, trim_a))
            os.symlink(os.path.abspath(b), os.path.join(tmp_dir, trim_b))
            file_list.append([(trim_a, trim_b)])
            basename_list.append("_".join(get_basename_noext(a).split("_")[:-1]))

    else:
        logging.info("The input dataset is considered as a single-end dataset.")
        for file in sorted(fastq_files):
            file_name = os.path.join(tmp_dir, get_basename_noext(file)) + ".{}".format(correct_ext(os.path.basename(file).split('.')[1:]))
            print(file_name)
            os.symlink(os.path.abspath(file), file_name)
            file_list.append(os.path.basename(file))
            basename_list.append(get_basename_noext(file))

    ret = dict()
    ret["paired"] = is_paired
    ret["inpath"] = tmp_dir
    
    ret["num_fastq"] = int(len(fastq_files) / (is_paired*2)) if is_paired else len(fastq_files)
    ret["file_list"] = file_list
    ret['basename_list'] = basename_list
    return ret

def run_NAMEseq_pipeline(param):
    import snakemake
    snakefile = os.path.join(os.path.dirname(__file__), "scripts/Snakefile.paired.py" if param["paired"] else "scripts/Snakefile.single.py")
    snakemake.snakemake(
        snakefile=snakefile,
        config={
            "input_path": param["inpath"],
            "index": param["--index"],
            "output_path": param["--outpath"],
            "num_threads" : param["--num_threads"],
            "basename_list": param["basename_list"],
            "fasta": param["--input_fasta"],
        },
        quiet=True
    )
    
def build_bowtie_index(in_fasta, ref_id):
    AT_only_fasta = in_fasta.split('.')[0]+'_AT_only.fasta'
    print(AT_only_fasta)
    logging.info("Building Bowtie Index")
    if not os.path.exists(in_fasta):
        logging.error("Input file does not exist.")
        sys.exit(1)
    ext = get_ext(in_fasta)
    if not ext.lower() in [ "fa", "fasta" ]:
        logging.error("Input file is not a FASTA file.")
        sys.exit(1)
    ref_path = os.path.join(os.path.dirname(__file__), "index")
    out_file_basename = os.path.join(ref_path, ref_id)

    if not os.path.exists(ref_path):
        os.mkdir(ref_path)

    cmd = "python scripts/fasta_to_AT_only.py {} {} ".format(in_fasta, AT_only_fasta)
    os.system(cmd)
    cmd = "bowtie2-build {} {} ".format(AT_only_fasta, out_file_basename)
    os.system(cmd)
    logging.info("Building '{}' index was finished!".format(ref_id))
    

def run(args):
    if args['quant']:
        print(os.path.dirname(__file__))
        if args['--num_threads'] is None:
            args['--num_threads'] = 4
        if args['--outpath'] is None:
            args['--outpath'] = os.path.join(os.getcwd(), "NAMEseq_output/")
        if args['--index'] is None or args['--index'] == "hs":
            args['--index'] = os.path.join(os.path.dirname(__file__), "index/hs")
        elif os.path.exists(os.path.join(os.path.dirname(__file__), "index", args['--index']+'.1.bt2')):
            args['--index'] = os.path.join(os.path.dirname(__file__), "index", args['--index'])
        else:
            logging.error("Reference file is not found!")
            sys.exit(1)

        if args['--index'].startswith("./"):
           args['--index'] = args['--index'][2:]
           
        logging.info("Starting quantification mode")
        logging.info("Collecting FASTQ files...")
        param = {**args, **collect_FASTQ_files(args['FILE'])}
        logging.info("Collected {} FASTQ files.".format(param["num_fastq"]))
        logging.info("Quantification has been finished.")
        logging.info("Running Snakemake pipline")

        run_NAMEseq_pipeline(param)

    if args['index']:
        build_bowtie_index(args['--input_fasta'], args['--index_name'])

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    args = docopt(__doc__, version='NAME-seq v0.1')
    run(args)