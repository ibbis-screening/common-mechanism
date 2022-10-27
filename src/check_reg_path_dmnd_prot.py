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
vax_ids = pd.read_csv(os.environ['DB_PATH'] + '/benign/vax_taxids', header=None)

if checkfile(file) != 1:
    exit(1)
diamond = readdmnd(file)
diamond = taxdist(diamond, reg_ids, vax_ids, query)

# trim down to the top hit for each region, ingnoring any top hits that are synthetic constructs
diamond2 = trimblast(diamond[diamond['subject tax ids']!="32630"])

reg_bac = 0
reg_vir = 0
reg_fung = 0

if diamond2['regulated'].sum(): # if ANY of the hits are regulated
    print("Regulated pathogen proteins: PRESENT")
    for gene in set(diamond2['subject acc.'][diamond2['regulated'] == True]):
        # go back to blast - the full set of hits
        # if it's a viral protein
        subset = diamond[(diamond['subject acc.'] == gene)]
        subset = subset.reset_index(drop=True)
        if subset['regulated'][0] == True:
            if subset['superkingdom'][0] == "Viruses":
                print("Regulated virus top hit: FLAG")
                reg_vir = 1
            elif subset['superkingdom'][0] == "Bacteria": # if top hit that isn't a synthetic construct is regulated
                print("Regulated bacteria top hit: FLAG")
                reg_bac = 1
            elif subset['superkingdom'][0] == "Fungi":
                print("Regulated fungi top hit: FLAG")
                reg_fung = 1
        n_reg = diamond['regulated'][diamond['subject acc.'] == gene].sum()
        n_total = len(diamond['regulated'][diamond['subject acc.'] == gene])
        if (n_reg < n_total):
            print("Gene " + gene + " found in both regulated and nonregulated organisms")
            print("Species: " + " ".join(set(diamond['species'][diamond['subject acc.'] == gene])))
        elif (n_reg == n_total):
            print("Gene " + gene + " found in only regulated organisms")
            print("Species: " + ", ".join(set(diamond['species'][diamond['subject acc.'] == gene])) + " (taxid: " + " ".join(map(str, set(diamond['subject tax ids'][diamond['subject acc.'] == gene]))) + ")")
        else:
            print("Gene: " + gene)
            print(diamond['regulated'][diamond['subject acc.'] == gene])
    hits = diamond2[diamond2['regulated']==True][['q. start', 'q. end']]  # print out the start and end coordinated on the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)


if reg_bac == 0:
    print("No regulated bacteria top hit: PASS")
if reg_vir == 0:
    print("No regulated virus top hit: PASS")
if reg_fung == 0:
    print("No regulated fungi top hit: PASS")
