#! /usr/bin/env python

#################################################################################
#check_benigh.py checks the output from hmmscan and prints to screen the results
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
#################################################################################
#Usage: 
#   python check_benign.py -i INPUT -s SEQUENCE -d DATABASE FOLDER 
#       -i, --input = input sample name (will check for sample.benign.hmmsearch file)
#       -s, --sequence = input sequence file
#       -d, --database = database folder location/path (will check for benign_annotations.csv) 
#################################################################################
from utils import *
import os, sys, argparse 
import pandas as pd
from itertools import chain
#import SeqIO

def check_for_benign(query, coords):
    # for each set of hits, need to pull out the coordinates covered by benign entries
    # overall >90% of regulated pathogen sub-sequences must be covered with benign content
    hmmsearch = query + ".benign.hmmsearch"
    qlen = 0
    if check_blastfile(hmmsearch) == 2:
        sys.stderr.write("\t...no housekeeping genes found")
        exit(0)
        
    hmmer = readhmmer(hmmsearch)
    for region in range(0, coords.shape[0]): # for each regulated pathogen region
        # look at only the hmmer hits that overlap with it
        htrim = hmmer[~((hmmer['ali from'] > coords['q. end'][region]) & (hmmer['ali to'] > coords['q. end'][region])) & ~((hmmer['ali from'] < coords['q. start'][region]) & (hmmer['ali to'] < coords['q. start'][region]))]
        if htrim.shape[0] > 0:
            htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'])
            if any(htrim['coverage'] > 0.90):
                htrim = htrim[htrim['coverage'] > 0.90]
                htrim = htrim.reset_index(drop=True)
                descriptions = []
                for row in range(htrim.shape[0]):
                    hit = htrim['target name'][row]
                    hit = hit.replace(".faa.final_tree.fa", "")
                    hit = hit.replace(".faa.final_tree.used_alg.fa", "")
                    descriptions.append(benign_desc['Annotation'][benign_desc['Model'] == hit])
                annot_string = "\n".join(str(v) for v in descriptions)
                sys.stderr.write("\t...housekeeping genes - >90% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved = PASS\n")
                sys.stderr.write(annot_string + "\n")
            else:
                sys.stderr.write("\t...housekeeping genes - <90% coverage achieved = FAIL\n")
        else:
            sys.stderr.write("\t...regulated region failed to clear\n")

    # annotate (and clear?) synbio parts
    blast = sys.argv[1] + ".benign.blastn"
    if check_blastfile(blast) == 2:
        sys.stderr.write("\t...No synbio parts hits\n")
    else:
        blastn = readblast(blast)
        for region in range(0, coords.shape[0]): # for each regulated pathogen region
            # look at only the hits that overlap with it
            htrim = blastn[~((blastn['q. start'] > coords['q. end'][region]) & (blastn['q. end'] > coords['q. end'][region])) & ~((blastn['q. start'] < coords['q. start'][region]) & (blastn['q. end'] < coords['q. start'][region]))]
            htrim = htrim.assign(coverage = abs(htrim['q. end'] - htrim['q. start']) / htrim['query length'])
            if any(htrim['coverage'] > 0.90):
                sys.stderr.write("\t...Synbio parts - >90% coverage achieved = PASS\n")
            else:
                sys.stderr.write("\t...Synbio parts - <90% coverage achieved = FAIL\n")

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="sample_name",
        required=True, help="Sample name")
    parser.add_argument("-s","--sequence", dest="seq_file",
        required=True, help="FASTA sequence file")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="Benign HMM database folder (must contain benign_annotations.csv)")
    args=parser.parse_args()

    #check input files
    if (not os.path.exists(args.seq_file)):
        sys.stderr.write("\t...sequence file does not exist\n")
        exit(1)
    if (not os.path.exists(args.db + "/benign_annotations.csv")):
        sys.stderr.write("\t...benign_annotations.csv does not exist\n")
        exit(1) 
    
    #Read in database file
    pd.set_option('max_colwidth',200)
    benign_desc = pd.read_csv(args.db + "/benign_annotations.csv")
    
    #Check for file - if exists, check for benign 
    if os.path.exists(args.sample_name + ".reg_path_coords.csv"):
        coords = pd.read_csv(args.sample_name + ".reg_path_coords.csv")
        check_for_benign(args.sample_name, coords)
    else:
        sys.stdout.write("\t...no regulated protein regions to clear\n")

    if os.path.exists(args.sample_name + ".reg_path_coords_nt.csv"):
        coords = pd.read_csv(args.sample_name + ".reg_path_coords_nt.csv")
        check_for_benign(args.sample_name, coords)
    else:
        sys.stdout.write("\t...no regulated nucleotide regions to clear\n")
    
if __name__ == "__main__":
    main()
