#! /usr/bin/env python

##################################################################################################
#check_reg_path_dmnd_prot.py checks results for regulated pathogen and prints any 
#matched coordinates
#   This file will ignore any synthetic constructs
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
##################################################################################################
#Usage:
#   python check_reg_path_dmnd_prot.py -i INPUT
#       -i, --input 
##################################################################################################
from utils import *
import os, sys, argparse
import pandas as pd
import taxoniq

def main(): 
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i","--input",dest="in_file", 
        required=True, help="Input query file (QUERY.nr.dmnd)")
    parser.add_argument("--benign-db", dest="benign_db",
        required=True,help="Benign HMM database folder (must contain vax_taxids file)")
    parser.add_argument("--biorisk-db", dest="biorisk_db",
        required=True,help="Biorisk HMM database folder (must contain reg_taxids file)")
    args=parser.parse_args() 
    
    #check input files
    if (not os.path.exists(args.in_file)):
        sys.stderr.write("\t...input query file %s does not exist\n" % args.in_file)
        exit(1)
    if (not os.path.exists(args.benign_db + "/vax_taxids")):
        sys.stderr.write("\t...benign db file %s does not exist\n" % (args.benign_db + "/vax_taxids"))
        exit(1)
    if (not os.path.exists(args.biorisk_db + "/reg_taxids")):
        sys.stderr.write("\t...biorisk db file %s does not exist\n" % (args.biorisk_db + "/reg_taxids"))
        exit(1)
    # read in files
    reg_ids = pd.read_csv(args.biorisk_db + "/reg_taxids", header=None)
    vax_ids = pd.read_csv(args.benign_db + "/vax_taxids", header=None)

    if check_blastfile(args.in_file) == 0:
        sys.stdout.write("\tERROR: Homology search has failed\n")
        exit(1)
    if check_blastfile(args.in_file) == 2:
        sys.stdout.write("\t...no protein hits\n")
        exit(1)
    blast = readblast(args.in_file)                  #function in utils.py
    blast = taxdist(blast, reg_ids, vax_ids) #function in utils.py

    # trim down to the top hit for each region, ignoring any top hits that are synthetic constructs
    blast2 = trimblast(blast[blast['subject tax ids']!="32630"])

    reg_bac = 0
    reg_vir = 0
    reg_fung = 0

    if blast2['regulated'].sum(): # if ANY of the trimmed hits are regulated
        sys.stdout.write("\t...regulated pathogen sequence: PRESENT")
        # print(blast[['subject acc.', 'regulated', 'genus', 'species']])
        # for each hit (subject acc) linked with at least one regulated taxid
        for gene in set(blast2['subject acc.'][blast2['regulated'] == True]): 
            # go back to blast - the full set of hits
            subset = blast[(blast['subject acc.'] == gene)]
            subset = subset.reset_index(drop=True)
            org = ""
            # if the top hit is found in regulated pathogens
            if subset['regulated'][0] == True:
                n_reg = blast['regulated'][blast['subject acc.'] == gene].sum()
                n_total = len(blast['regulated'][blast['subject acc.'] == gene])
                # if some of the organisms with this gene aren't regulated, say so
                if (n_reg < n_total):
                    sys.stdout.write("\t\t --> %s found in both regulated and non-regulated organisms\n" % gene)
                    sys.stdout.write("\t...species: %s\n" % (" ".join(set(blast['species'][blast['subject acc.'] == gene])))) 
                    # could explicitly list which are and aren't regulated?
                    # otherwise, raise a flag and say which superkingdom the flag belongs to
                elif (n_reg == n_total):
                    if subset['superkingdom'][0] == "Viruses":
                        reg_vir = 1
                        org = "virus"
                    elif subset['superkingdom'][0] == "Bacteria": 
                        reg_bac = 1
                        org = "bacteria"
                    elif subset['kingdom'][0] == "Fungi":
                        org = "fungi"
                        reg_fung = 1
                    sys.stdout.write("\t...%s\n" % (subset['superkingdom'][0]))
                    sys.stdout.write("\t\t --> %s found in only regulated organisms: FLAG (%s)\n" % (gene, org))
                    sys.stdout.write("\t...species: %s (taxid: %s)\n" % ((", ".join(set(blast['species'][blast['subject acc.'] == gene]))), (" ".join(map(str, set(blast['subject tax ids'][blast['subject acc.'] == gene]))))))
                else: # something is wrong, n_reg > n_total
                    sys.stdout.write("\t...gene: %s\n" % gene)
                    sys.stdout.write("%s\n" % (blast['regulated'][blast['subject acc.'] == gene]))
        hits = blast2[blast2['regulated']==True][['q. start', 'q. end']]  # print out the start and end coordinates of the query sequence
        #Create output file 
        sample_name = args.in_file[0:-8]
        hits.to_csv(sample_name + ".reg_path_coords.csv", index=False)

    if reg_vir == 0 and reg_bac == 0 and reg_fung == 0:
        sys.stdout.write("\t\t --> no regulated pathogen top hit: PASS\n")

if __name__ == "__main__":
    main()
