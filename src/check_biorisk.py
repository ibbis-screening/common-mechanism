#! /usr/bin/env python3

##############################################################################
#check_biorisk.py checks the output from hmmscan and prints to screen the results
#
#Copyright (C) 2022-2023 NTI|Bio 
#This file is part of the CommonMechanism 
##############################################################################
# Usage:
#  python check_biorisk.py -i INPUT.biorisk.hmmscan -d databases/biorisk_db/ 
##############################################################################
from utils import *
import os, sys, argparse 
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="in_file",
        required=True, help="Input file - hmmscan output file")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="HMM folder (must contain biorisk_annotations.csv)")
    args = parser.parse_args()
    
    #check input files
    if (not os.path.exists(args.in_file)):
        sys.stderr.write("\t...input file does not exist\n") 
        exit(1) 
    if (not os.path.exists(args.db + "/biorisk_annotations.csv")):
        sys.stderr.write("\t...biorisk_annotations.csv does not exist\n")
        exit(1)
    
    #Specify input file and read in database file 
    in_file = args.in_file
    # sys.stdout.write("\t...checking %s\n" % in_file) 

    lookup = pd.read_csv(args.db + "/biorisk_annotations.csv")
    lookup.fillna(False, inplace=True)
    # print(lookup)

    # read in HMMER output and check for valid hits
    res = check_blastfile(in_file)
    if res == 0:
        sys.stdout.write("\t...ERROR: biorisk search results empty\n")
    if res == 1:
        hmmer = readhmmer(in_file)
        keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-20]
        hmmer = hmmer.iloc[keep1,:]
        hmmer = trimhmmer(hmmer)
        hmmer['description'] = ''
        hmmer['Must flag'] = False
        hmmer = hmmer.reset_index(drop=True)
        # print(hmmer)
        for model in range(hmmer.shape[0]):
            name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
            # print(name_index)
            hmmer.loc[model, 'description'] = lookup.iloc[name_index[0], 1]
            hmmer.loc[model, 'Must flag'] = lookup.iloc[name_index[0], 2]
        if hmmer.shape[0] > 0:
            if (sum(hmmer['Must flag']) > 0):
                for region in hmmer.index[hmmer['Must flag'] != 0]:
                    if (hmmer['ali from'][region] > hmmer['qlen'][region]):
                        hmmer['ali from'][region] = divmod(hmmer['ali from'][region], hmmer['qlen'][region])[0]
                        hmmer['ali to'][region] = divmod(hmmer['ali to'][region], hmmer['qlen'][region])[0]
                    sys.stdout.write("\t\t --> Biorisks: Regulated gene in bases " + str(hmmer['ali from'][region]) + " to " + str(hmmer['ali to'][region]) + ": FLAG\n\t\t     Gene: " + ", ".join(set(hmmer['description'][hmmer['Must flag'] == True])) + "\n")
            else:
                sys.stdout.write("\t\t --> Biorisks: Regulated genes not found, PASS\n")
            if (sum(hmmer['Must flag']) != hmmer.shape[0]):
                for region in hmmer.index[hmmer['Must flag'] == 0]:
                    sys.stdout.write("\t\t --> Virulence factor found in bases " + str(hmmer['ali from'][region]) + " to " + str(hmmer['ali to'][region]) + ", WARNING\n\t\t     Gene: " + ", ".join(set(hmmer['description'][hmmer['Must flag'] == False])) + "\n")
        else: 
            sys.stdout.write("\t\t --> Biorisks: no significant hits detected, PASS\n")
    if res == 2:
        sys.stdout.write("\t\t --> Biorisks: no hits detected, PASS\n")

if __name__ == "__main__":
    main()
