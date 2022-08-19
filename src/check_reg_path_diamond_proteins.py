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
file = query + ".nr.diamond"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

diamond = readdmnd(file)
diamond = taxdist(diamond, reg_ids, query)
#print(diamond.iloc[:,0:6])
diamond = trimblast(diamond)
#print(diamond)

# ignore synthetic constructs when deciding whether to flag
#print(diamond)

if diamond['regulated'].sum(): # if ANY of the hits are regulated
    print("Regulated pathogen proteins: PRESENT")
    if "Viruses" in set(diamond['superkingdom'][diamond['regulated'] == True]):
        print("Regulated virus: FLAG")
    elif diamond['regulated'][diamond['subject tax ids']!="32630"][0] == True: # if top hit that isn't a synthetic construct is regulated
        print("Regulated bacteria top hit: FLAG")
#    hits = diamond[diamond['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
    hits = diamond[diamond['regulated']==True]  # print out the start and end coordinated on the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv", index=False)
else:
	print("Regulated pathogen proteins: PASS")












# specify which parts of the order aren't covered by a protein



#
#
## read in nt screening
#
#file = sys.argv[1] + ".nt.blastn"
#
#blast = taxdist(file, reg_ids, query)
##print(blast)
#blast = trimblast(blast)
##print(blast)
#
#if blast['regulated'].sum():
#    print("Regulated pathogens (nt): FLAG")
#    if "Viruses" in set(blast['superkingdom'][blast['regulated'] == True]):
#        print("Regulated virus (nt)")
#    elif blast['regulated'][blast['subject tax ids']!="32630"][0] == True: # if top hit that isn't a synthetic construct is regulated
#        print("Regulated bacteria top hit: FLAG")
#    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
##    print(hits)
#    hits.to_csv(sys.argv[1] + ".reg_path_coords_nt.csv", index=False)
#else:
#    print("Regulated pathogens (nt): PASS")
#
## screen any nt regions of concern
#bedtools getfasta -fi query.fasta -bed test.bed -fo test.fa.out
#
#$ cat test.bed
#chr1 5 10
