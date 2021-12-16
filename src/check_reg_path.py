from utils import *
import os, sys
import pandas as pd
import taxoniq

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

file = sys.argv[1] + ".nr.blastx"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

blast = taxdist(file, reg_ids)
#print(blast)
blast = trimblast(blast)
#print(blast)

if blast['regulated'].sum():
#    print("Regulated pathogens: FLAG")
    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
#    print(hits)
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv")
else:
	print("Regulated pathogens: PASS")
