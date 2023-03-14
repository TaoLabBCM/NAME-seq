#!/usr/bin/python

import sys
import re

fastq_file = sys.argv[1]
converted_sample = sys.argv[2]
with open(fastq_file, "r") as f_original:
    with open(converted_sample, "w") as f_converted:
        n = 0
        for line in f_original:
            n += 1
            if n%4 == 2:
                line = re.sub('G', 'A', line)
                line = re.sub('C', 'T', line)
            f_converted.write(line)