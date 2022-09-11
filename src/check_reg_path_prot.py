from utils import *
import os, sys
import pandas as pd
import taxoniq

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates
# we have opted to ignore any "synthetic construct" hits
# usage: python -m check_reg_path ${name}

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

query = sys.argv[1]

# read in protein screening
file = query + ".nr.blastx"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

blast = readblast(file)
blast = taxdist(blast, reg_ids, query)
#print(blast)
#print(blast['regulated'])

# trim down to the top hit for each region, ingnoring any top hits that are synthetic constructs
blast2 = trimblast(blast[blast['subject tax ids']!="32630"])
#print(blast2)
#print(blast.iloc[:,10:17])

# ignore synthetic constructs when deciding whether to flag

reg_bac = 0
reg_vir = 0

if blast2['regulated'].sum(): # if ANY of the trimmed hits are regulated
#    print(blast2)
    print("Regulated pathogen proteins: PRESENT")
    for gene in set(blast2['subject acc.'][blast2['regulated'] == True]): # for each gene with at least one regulated hit
        # go back to blast - the full set of hits
        # if it's a viral protein
        subset = blast[(blast['subject tax ids'] != "32630") & (blast['subject acc.'] == gene)]
        subset = subset.reset_index(drop=True)
        if "Viruses" in set(subset['superkingdom']):
            # if the top hit is both viral and regulated
            if subset['superkingdom'][0] == "Viruses":
                if subset['regulated'][0] == True:
                    print("Regulated virus top hit: FLAG")
                    reg_vir = 1
            else:
                print("Regulated virus in hits: COND FLAG")
        elif subset['regulated'][0] == True: # if top hit that isn't a synthetic construct is regulated
            print("Regulated bacteria top hit: FLAG")
            reg_bac = 1
        n_reg = blast['regulated'][blast['subject acc.'] == gene].sum()
        n_total = len(blast['regulated'][blast['subject acc.'] == gene])
        if (n_reg < n_total):
            print("Gene " + gene + " found in both regulated and nonregulated organisms")
            print("Species: " + " ".join(set(blast['species'][blast['subject acc.'] == gene])))
        elif (n_reg == n_total):
            print("Gene " + gene + " found in only regulated organisms")
            print("Species: " + " ".join(set(blast['species'][blast['subject acc.'] == gene])) + " (taxid: " + " ".join(map(str, set(blast['subject tax ids'][blast['subject acc.'] == gene]))) + ")")
        else:
            print("Gene: " + gene)
            print(blast['regulated'][blast['subject acc.'] == gene])
    hits = blast2[blast2['regulated']==True][['q. start', 'q. end']]  # print out the start and end coordinates of the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)

if reg_bac == 0:
    print("No regulated bacteria top hit: PASS")
if reg_vir == 0:
    print("No regulated virus top hit: PASS")


