#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks whether there are any hits to nr for a query. If there aren't any over a given
significance level, prints the whole sequence to a noncoding query file. If there are hits, fetches
the nucleotide regions between these hits and singles them out for nucleotide screening

Usage:
    fetch_nc_bits.py query_name fasta_file_path
"""
import sys, shutil
import re
from Bio import SeqIO
import argparse
#from commec.utils import *
from commec.tools.blast_tools import BlastTools
from commec.utils.file_utils import FileTools

    #query = sys.argv[1]
    #f_file = sys.argv[2]

def fetch_noncoding_regions(nr_output_file : str, cleaned_fasta_file_path : str):
    """
    Takes the output of a blastx result from NR, and the query fasta file, and fetches a list of 
    """
    nc_bits = 0
    query = nr_output_file
    f_file = cleaned_fasta_file_path

    # check if the nr hits file is empty
    if FileTools.is_empty(query):
        nc_bits = "all"
    elif not FileTools.has_hits(query):
        sys.stdout.write("\t...no hits to the nr database\n")
        nc_bits = "all"
    # if not, check whether any of the hits has an E-value > 1e-30
    # if so, find the start and end of those hits and use these to get the coordinates of non-coding regions
    else:
        blast = BlastTools.readblast(query)
        blast = BlastTools.trimblast(blast)
        blast = blast[blast['% identity'] >= 90]
        if blast.shape[0] > 0:
        # find noncoding bits
            sys.stdout.write("\t...protein hits found, fetching nt regions not covered by a 90% ID hit or better\n")
            hits = []
            for i in range(blast.shape[0]):
                pair = [blast['q. start'][i], blast['q. end'][i]]
                pair.sort()
                hits.append(pair)
            hits = sorted(hits, key=lambda x: x[0])

            nc_bits = []
            if hits[0][0] >50:
                nc_bits.append([1,hits[0][0]])
            for i in range(len(hits)-1):
                if hits[i][1] < (hits[i+1][0] - 49): # if there's a noncoding region of >=50 between hits
                    nc_bits.append([hits[i][1], hits[i+1][0]])
        else:
            sys.stdout.write("\t...protein hits all low percent identity (<90%) - screening entire sequence\n")
            nc_bits = "all"

    # fetch noncoding sequences

    outfile = re.sub(".nr.*", "", query) + '.noncoding.fasta'

    def fetch_sequences(seqid, nc_bits, f_file, outfile):
        tofetch = []
        for (start, stop) in nc_bits:
            tofetch.append((seqid, start, stop))

        if tofetch:
            with open(f_file, "r") as fasta_file:
                records = list(SeqIO.parse(fasta_file, "fasta"))
                sequences = []
                for (seqid, start, stop) in tofetch:
                    start = int(start)
                    stop = int(stop)
                    for record in records:
                        sequence = record.seq[start - 1 : stop]  # Adjust start to 0-based index
                        sequences.append(f">{seqid} {start}-{stop}\n{sequence}\n")
                        break
            with open(outfile, "w", encoding = "utf-8") as output_file:
                output_file.writelines(sequences)

    if nc_bits == "all":
        shutil.copyfile(f_file, outfile)
    elif nc_bits == []: # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
        sys.stdout.write("\t\t --> no noncoding regions >= 50 bases found, skipping nt scan\n")
    else:
        seqid = blast.iloc[0, 0]
        fetch_sequences(seqid, nc_bits, f_file, outfile)

def main():
    '''
    Wrapper for parsing arguments direction to fetch_nc_bits if called as main.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="in_file",
        required=True, help="Input query file path")
    parser.add_argument("-f","--fasta", dest="fasta_file",
        required=True, help="HMM folder (must contain biorisk_annotations.csv)")
    args = parser.parse_args()
    fetch_noncoding_regions(args.in_file, args.fasta_file_path)

if __name__ == "__main__":
    main()