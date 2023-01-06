from utils import *
import os, sys
import pandas as pd

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates
# we have opted to ignore any "synthetic construct" hits
# usage: python -m check_reg_path ${name}

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

file = sys.argv[1]

# read in protein screening
reg_ids = pd.read_csv(os.environ['DB_PATH'] + '/biorisk/reg_taxids', header=None)
vax_ids = pd.read_csv(os.environ['DB_PATH'] + '/benign/vax_taxids', header=None)

if check_blastfile(file) != 1:
    exit(1)
blast = readblast(file)
blast = taxdist(blast, reg_ids, vax_ids)

# trim down to the top hit for each region, ingnoring any top hits that are synthetic constructs
blast2 = trimblast(blast[blast['subject tax ids']!="32630"])

reg_bac = 0
reg_vir = 0
reg_fung = 0

if blast2['regulated'].sum(): # if ANY of the trimmed hits are regulated
#    print(blast2)
    print("Regulated pathogen proteins: PRESENT")
    print(blast[['subject acc.', 'regulated', 'genus', 'species']])
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
                print("Gene " + gene + " found in both regulated and non-regulated organisms")
                print("Species: " + " ".join(set(blast['species'][blast['subject acc.'] == gene]))) # could explicitly list which are and aren't regulated?
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
                print(subset['superkingdom'][0])
                print("Gene " + gene + " found in only regulated organisms: FLAG (" + org + ")")
                print("Species: " + ", ".join(set(blast['species'][blast['subject acc.'] == gene])) + " (taxid: " + " ".join(map(str, set(blast['subject tax ids'][blast['subject acc.'] == gene]))) + ")")
            else: # something is wrong, n_reg > n_total
                print("Gene: " + gene)
                print(blast['regulated'][blast['subject acc.'] == gene])
    hits = blast2[blast2['regulated']==True][['q. start', 'q. end']]  # print out the start and end coordinates of the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)

if reg_vir == 0 and reg_bac == 0 and reg_fung == 0:
    print("No regulated pathogen top hit: PASS")


