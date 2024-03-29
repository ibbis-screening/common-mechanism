#!/usr/bin/env python
#Jennifer Lu, jlu26@jhmi.edu
#2/3/2022
#
#This program splits a multi-fasta file based on number of sequences per file
#Output will be input.#.fa

import argparse
import os
import sys
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="i_file",
        required=True, help="multi-FASTA file to split")
    parser.add_argument("-n","--num", dest="num_seqs",
        type=int,
        required=True, help="Number of sequences per file (min)")
    args = parser.parse_args()

    basename = os.path.splitext(args.i_file)[0]
    basename = os.path.basename(basename)
    count_curr = 0
    count_total = 0
    num_splits = 0
    sys.stdout.write(f"\t{count_total} sequences printed ({num_splits} splits)")
    sys.stdout.flush()
    for record in SeqIO.parse(args.i_file,"fasta"):
        if count_curr == 0:
            num_splits += 1
            o_file = open(basename + "." + str(num_splits) + ".fa" , 'w')
        SeqIO.write(record, o_file, "fasta")
        count_curr += 1
        count_total += 1
        if count_total % 10000 == 0:
            sys.stdout.write("\r\t{count_total} sequences printed ({num_splits} splits)")
            sys.stdout.flush()
        if count_curr == args.num_seqs:
            count_curr = 0
            o_file.close()

    sys.stdout.write(f"\t{count_total} sequences printed ({num_splits} splits)\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
