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
#print(blast['regulated'])

# trim down to the top hit for each region
blast2 = trimblast(blast)
#print(blast.iloc[:,10:17])

# ignore synthetic constructs when deciding whether to flag

if blast2['regulated'].sum(): # if ANY of the hits are regulated
    print("Regulated pathogen proteins: PRESENT")
    for gene in set(blast2['subject acc.'][blast2['regulated'] == True]):
        if "Viruses" in set(blast['superkingdom'][blast['subject acc.'] == gene]):
            print("Regulated virus: FLAG")
        elif blast['regulated'][blast['subject tax ids']!="32630"][0] == True: # if top hit that isn't a synthetic construct is regulated
            print("Regulated bacteria top hit: FLAG")
        n_reg = blast['regulated'][blast['subject acc.'] == gene].sum()
        n_total = len(blast['regulated'][blast['subject acc.'] == gene])
        if (n_reg < n_total):
            print("Gene " + gene + " found in both regulated and nonregulated organisms: COND FLAG")
            print(" ".join(blast['species'][blast['subject acc.'] == gene]))
        elif (n_reg == n_total):
            print("Gene " + gene + " found in only regulated organisms: FLAG")
            print(" ".join(blast['species'][blast['subject acc.'] == gene]))
        else:
            print("Gene: " + gene)
            print(blast['regulated'][blast['subject acc.'] == gene])
#    hits = diamond[diamond['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
    hits = blast[blast['regulated']==True]  # print out the start and end coordinated on the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)
else:
    print("Regulated pathogen proteins: PASS")



#if blast['regulated'].sum(): # if ANY of the hits are regulated
#    print("Regulated pathogens: PRESENT")
#    if "Viruses" in set(blast['superkingdom'][blast['regulated'] == True]):
#        print("Regulated virus: FLAG")
#    elif blast['regulated'][blast['subject tax ids']!="32630"][0] == True:
#        print("Regulated bacteria top hit: FLAG")
#    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
##    print(hits)
#    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)
#else:
#    print("Regulated pathogens: PASS")
