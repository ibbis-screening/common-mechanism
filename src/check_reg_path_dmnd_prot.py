from utils import *
import os, sys
import pandas as pd
import taxoniq

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates
# we have opted to ignore any "synthetic construct" hits
# usage: python src/check_reg_path_diamond.py ${file prefix}

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)
    
query = sys.argv[1]

# read in protein screening results
file = query + ".nr.dmnd"
reg_ids = pd.read_csv(os.environ['DB_PATH'] + '/biorisk/reg_taxids', header=None)

if checkfile(file) == 1:
    diamond = readdmnd(file)
    diamond = readdmnd(file)
    diamond = taxdist(diamond, reg_ids, query)
    diamond2 = trimblast(diamond)

    # ignore synthetic constructs when deciding whether to flag
    #print(diamond)

    if diamond['regulated'].sum(): # if ANY of the hits are regulated
        print("Regulated pathogen proteins: PRESENT")
        for gene in set(diamond2['subject acc.'][diamond2['regulated'] == True]):
            if "Viruses" in set(diamond['superkingdom'][diamond['subject acc.'] == gene]):
                print("Regulated virus: FLAG")
            elif diamond['regulated'][diamond['subject tax ids']!="32630"][0] == True: # if top hit that isn't a synthetic construct is regulated
                print("Regulated bacteria top hit: FLAG")
            n_reg = diamond['regulated'][diamond['subject acc.'] == gene].sum()
            n_total = len(diamond['regulated'][diamond['subject acc.'] == gene])
            if (n_reg < n_total):
                print("Gene " + gene + " found in both regulated and nonregulated organisms: COND FLAG")
            elif (n_reg == n_total):
                print("Gene " + gene + " found in only regulated organisms: FLAG")
            else:
                print("Gene: " + gene)
                print(diamond['regulated'][diamond['subject acc.'] == gene])
    #    hits = diamond[diamond['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
        hits = diamond[diamond['regulated']==True]  # print out the start and end coordinated on the query sequence
        hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)
    else:
        print("Regulated pathogen proteins: PASS")
elif checkfile(file) == 2:
    print("No protein hits found")
else:
    print("Unexpected result from DIAMOND nr screen")

