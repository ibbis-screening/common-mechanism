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

# read in protein screening
file = sys.argv[1] + ".nr.blastx"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

blast = readblast(file)
blast = taxdist(blast, reg_ids, sys.argv[1])
#print(blast)

# trim down to the top hit for each region
blast = trimblast(blast)
#print(blast)

# ignore synthetic constructs when deciding whether to flag

if blast['regulated'].sum(): # if ANY of the hits are regulated
    print("Regulated pathogens: PRESENT")
    if "Viruses" in set(blast['superkingdom'][blast['regulated'] == True]):
        print("Regulated virus: FLAG")
    elif blast['regulated'][blast['subject tax ids']!="32630"][0] == True:
        print("Regulated bacteria top hit: FLAG")
    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
#    print(hits)
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)
else:
	print("Regulated pathogens: PASS")

# read in nt screening

file = sys.argv[1] + ".nt.blastn"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

blast = taxdist(file, reg_ids, sys.argv[1])
#print(blast)
blast = trimblast(blast)
#print(blast)

if blast['regulated'].sum():
    print("Regulated pathogens (nt): FLAG")
    if "Viruses" in set(blast['superkingdom'][blast['regulated'] == True]):
        print("Regulated virus (nt)")
    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
#    print(hits)
    hits.to_csv(sys.argv[1] + ".reg_path_coords_nt.csv", index=False)
else:
    print("Regulated pathogens (nt): PASS")
