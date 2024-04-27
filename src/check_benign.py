#! /usr/bin/env python3

#################################################################################
#check_benigh.py checks the output from hmmscan and prints to screen the results
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
#################################################################################
#Usage: 
#   python check_benign.py -i INPUT -s SEQUENCE -d DATABASE FOLDER 
#       -i, --input = input sample name (will check for sample.benign.hmmscan file)
#       -s, --sequence = input sequence file
#       -d, --database = database folder location/path (will check for benign_annotations.tsv) 
#################################################################################
from utils import *
import os, sys, argparse 
import pandas as pd

def check_for_benign(query, coords, benign_desc):
        
    cleared = [0] * coords.shape[0]
    
    # PROTEIN HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    hmmscan = query + ".benign.hmmscan"
    if not has_hits(hmmscan):
        sys.stdout.write("\t...no housekeeping protein hits\n")
    else:
        hmmer = readhmmer(hmmscan)
        hmmer = hmmer[hmmer['E-value'] < 1e-20]
        # print(hmmer)
        for region in range(0, coords.shape[0]): # for each regulated pathogen region
            # look at only the hmmer hits that overlap with it
            htrim = hmmer[~((hmmer['ali from'] > coords['q. end'][region]) & (hmmer['ali to'] > coords['q. end'][region])) & ~((hmmer['ali from'] < coords['q. start'][region]) & (hmmer['ali to'] < coords['q. start'][region]))]
            if htrim.shape[0] > 0:
                # htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'])
                # if any(htrim['coverage'] > 0.80):
                #     # print(htrim)
                #     htrim = htrim[htrim['coverage'] > 0.80]
                #     htrim = htrim.reset_index(drop=True)
                #     descriptions = []
                #     for row in range(htrim.shape[0]):
                #         hit = htrim['target name'][row]
                #         descriptions.append(hit + ": " + str(*benign_desc['Description'][benign_desc['ID'] == hit]) + "\n")
                #     annot_string = "\n".join(str(v) for v in descriptions)
                #     sys.stdout.write("\t\t -->Housekeeping proteins - >80% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved = PASS\n")
                #     sys.stdout.write("\t\t   " + annot_string)
                #     cleared[region] = 1
                htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']))
                if any(htrim['qlen'] - htrim['coverage'] < 50):
                    # print(htrim)
                    htrim = htrim[htrim['coverage'] > 0.80]
                    htrim = htrim.reset_index(drop=True)
                    descriptions = []
                    # for row in range(htrim.shape[0]):
                    for row in [0]: # just print the top hit
                        hit = htrim['target name'][row]
                        # print(benign_desc['Description'][benign_desc['ID'] == hit])
                        descriptions.append(hit + ": " + str(*benign_desc['Description'][benign_desc['ID'] == hit]) + " (E-value: " + '{:.3g}'.format(htrim['E-value'][row], 3) + ")\n")
                    annot_string = "\n".join(str(v) for v in descriptions)
                    sys.stdout.write("\t\t -->Housekeeping proteins covering " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " = PASS\n")
                    sys.stdout.write("\t\t   " + annot_string)
                    cleared[region] = 1
                else:
                    # print(htrim)
                    sys.stdout.write("\t\t -->Housekeeping proteins - not enough coverage = FAIL\n")
                
    # RNA HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    cmscan = query + ".benign.cmscan"
    if not has_hits(cmscan):
        sys.stdout.write("\t...no benign RNA hits\n")
    else:
        cmscan = readcmscan(cmscan)
        # print(cmscan)
        for region in range(0, coords.shape[0]): # for each regulated pathogen region
            # look at only the cmscan hits that overlap with it
            qlen = abs(coords['q. start'][region] - coords['q. end'][region])
            # filter hits for ones that overlap with the regulated region
            htrim = cmscan[~((cmscan['seq from'] < coords['q. start'][region]) & (cmscan['seq to'] < coords['q. start'][region])) & ~((cmscan['seq from'] > coords['q. end'][region]) & (cmscan['seq to'] > coords['q. end'][region]))]
            # print(htrim)
            if htrim.shape[0] > 0:
                # percent coverage based method
                # htrim = htrim.assign(coverage = abs(htrim['seq to'] - htrim['seq from']) / qlen)
                # if any(htrim['coverage'] > 0.80):
                #     htrim = htrim[htrim['coverage'] > 0.80]
                #     htrim = htrim.reset_index(drop=True)
                #     descriptions = []
                #     for row in range(htrim.shape[0]):
                #         hit = htrim['target name'][row]
                #         descriptions.append(hit)
                #     annot_string = "\n\t...".join(str(v) for v in descriptions)
                #     sys.stdout.write("\t\t -->Housekeeping RNAs - >80% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved: PASS\n")
                #     sys.stdout.write("\t\t   RNA family: " + annot_string + "\n")
                #     cleared[region] = 1
                # else:
                #     sys.stdout.write("\t\t -->Housekeeping RNAs - <80% coverage achieved = FAIL\n")
                
                # bases unaccounted for based method
                htrim = htrim.assign(coverage = qlen - abs(htrim['seq to'] - htrim['seq from']))
                if any(htrim['coverage'] < 50):
                    htrim = htrim[htrim['coverage'] < 50]
                    htrim = htrim.reset_index(drop=True)
                    descriptions = []
                    for row in range(htrim.shape[0]):
                        hit = htrim['target name'][row]
                        descriptions.append(hit)
                    annot_string = "\n\t...".join(str(v) for v in descriptions)
                    sys.stdout.write("\t\t -->Housekeeping RNAs - <50 bases unaccounted for: PASS\n")
                    sys.stdout.write("\t\t   RNA family: " + annot_string + "\n")
                    cleared[region] = 1
                else:
                    sys.stdout.write("\t\t -->Housekeeping RNAs - >50 bases unaccounted for = FAIL\n")

    # SYNBIO HITS
    # annotate and clear benign nucleotide sequences
    blast = query + ".benign.blastn"
    if not has_hits(blast):
        sys.stdout.write("\t...no Synbio sequence hits\n")
    else:
        blastn = readblast(blast) # synbio parts
        blastn = trimblast(blastn)
        blastn = tophits(blastn)
        # print(blastn)
        for region in range(0, coords.shape[0]): # for each regulated pathogen region
            htrim = blastn[~((blastn['q. start'] > coords['q. end'][region]) & (blastn['q. end'] > coords['q. end'][region])) & ~((blastn['q. start'] < coords['q. start'][region]) & (blastn['q. end'] < coords['q. start'][region]))]
            if any(htrim['q. coverage'] > 0.80):
                htrim = htrim[htrim['q. coverage'] > 0.80]
                htrim = htrim.reset_index(drop=True)
                descriptions = []
                for row in range(htrim.shape[0]):
                    hit = htrim['subject title'][row]
                    descriptions.append(hit)
                annot_string = "\n\t\t   ".join(str(v) for v in descriptions)
                sys.stdout.write("\t\t -->Synbio sequences - >80% coverage achieved = PASS\n")
                sys.stdout.write("\t\t   Synbio parts: " + annot_string + "\n")
                cleared[region] = 1
            else:
                sys.stdout.write("\t\t -->Synbio sequences - <80% coverage achieved = FAIL\n")
            
    for region in range(0, coords.shape[0]):
        if cleared[region] == 0:
            sys.stdout.write("\t\t -->Regulated region at bases " + str(int(coords.iloc[region, 0])) + " to "  + str(int(coords.iloc[region, 1])) + " failed to clear: FLAG\n")
    if sum(cleared) == len(cleared):
        sys.stdout.write("\n\t\t -->all regulated regions cleared: PASS\n")

def main(): 
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input", dest="sample_name",
        required=True, help="Sample name")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="Benign HMM database folder (must contain benign_annotations.tsv)")
    args=parser.parse_args()

    if (not os.path.exists(args.db + "/benign_annotations.tsv")):
        sys.stderr.write("\t...benign_annotations.tsv does not exist\n")
        exit(1) 
    
    #Read in database file
    pd.set_option('max_colwidth',200)
    benign_desc = pd.read_csv(args.db + "/benign_annotations.tsv", sep="\t")
        
    #Check for file - if exists, check for benign 
    if os.path.exists(args.sample_name + ".reg_path_coords.csv"):
        coords = pd.read_csv(args.sample_name + ".reg_path_coords.csv")
        if coords.shape[0] == 0:
            sys.stdout.write("\t...no regulated regions to clear\n")
            exit(0)
        coords.sort_values(by=['q. start'], inplace=True)
        coords.reset_index(drop=True, inplace=True)
        check_for_benign(args.sample_name, coords, benign_desc)
    else:
        sys.stdout.write("\t...no regulated regions to clear\n")
    
if __name__ == "__main__":
    main()
