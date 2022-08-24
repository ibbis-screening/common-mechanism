from utils import *
import os, sys
import pandas as pd
import taxoniq

cutoff=200 # min length of a region not covered by a protien for follow-up nt scanning

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates
# we have opted to ignore any "synthetic construct" hits

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

# read in nucleotide screening

query = sys.argv[1]
file = $query + ".nt.blastn"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

blast = taxdist(file, reg_ids, $query)
#print(blast)
blast = trimblast(blast)
#print(blast)

if blast['regulated'].sum():
    print("Regulated pathogens (nt): FLAG")
    if "Viruses" in set(blast['superkingdom'][blast['regulated'] == True]):
        print("Regulated virus (nt)")
    if "Bacteria" in set(blast['superkingdom'][blast['regulated'] == True]):
        print("Regulated bacteria (nt)")
    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinates of the query sequence
#    print(hits)
    hits.to_csv(sys.argv[1] + ".reg_path_coords_nt.csv", index=False)
else:
    print("Regulated pathogens (nt): PASS")
