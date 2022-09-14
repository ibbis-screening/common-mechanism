from utils import *
import os, sys
import pandas as pd
from itertools import chain
#import SeqIO

#Input parameter error checking 
if len(sys.argv) < 2:
    sys.stdout.write("\tERROR: Please provide two files: one sequence name, one FASTA for the sequence\n")
    exit(1)
    
def check_for_benign(query, coords):
        # for each set of hits, need to pull out the coordinates covered by benign entries
        # overall >90% of regulated pathogen sub-sequences must be covered with benign content
        hmmsearch = sys.argv[1] + ".benign.hmmsearch"
        qlen = 0
        if checkfile(hmmsearch) == 2:
            print("No housekeeping genes found")
            exit(0)
        
        hmmer = readhmmer(hmmsearch)
#        print(hmmer)
        qlen = hmmer['qlen'][0]
        for region in range(0, coords.shape[0]): # for each regulated pathogen region
            # look at only the hmmer hits that overlap with it
            htrim = hmmer[~((hmmer['ali from'] > coords['q. end'][region]) & (hmmer['ali to'] > coords['q. end'][region])) & ~((hmmer['ali from'] < coords['q. start'][region]) & (hmmer['ali to'] < coords['q. start'][region]))]
            if htrim.shape[0] > 0:
                htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'])
    #            htrim['coverage'] = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'][0]
                if any(htrim['coverage'] > 0.90):
#                    print(" ".join(list(chain(*list(zip(htrim['target name'][htrim['coverage'] > 0.90], htrim['description of target'][htrim['coverage'] > 0.90]))))))
                    print("Housekeeping genes - >90% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved by " + " ".join(list(chain(*list(zip(htrim['target name'][htrim['coverage'] > 0.90], htrim['description of target'][htrim['coverage'] > 0.90]))))) + " = PASS")
                else:
                    print("Housekeeping genes - <90% coverage achieved = FAIL")
            else:
                print("Regulated region failed to clear")

        # annotate (and clear?) synbio parts
        blast = sys.argv[1] + ".benign.blastn"
        if checkfile(blast) == 2:
            print("No synbio parts hits")
        else:
            blastn = readblast(blast)
            for region in range(0, coords.shape[0]): # for each regulated pathogen region
                # look at only the hits that overlap with it
                htrim = blastn[~((blastn['q. start'] > coords['q. end'][region]) & (blastn['q. end'] > coords['q. end'][region])) & ~((blastn['q. start'] < coords['q. start'][region]) & (blastn['q. end'] < coords['q. start'][region]))]
                htrim = htrim.assign(coverage = abs(htrim['q. end'] - htrim['q. start']) / htrim['query length'])
                if any(htrim['coverage'] > 0.90):
                    print("Synbio parts - >90% coverage achieved = PASS")
                else:
                    print("Synbio parts - <90% coverage achieved = FAIL")
    


if os.path.exists(sys.argv[1] + ".reg_path_coords.csv"):
    coords = pd.read_csv(sys.argv[1] + ".reg_path_coords.csv")
    check_for_benign(sys.argv[1], coords)
else:
    print("No regulated protein regions to clear")

if os.path.exists(sys.argv[1] + ".reg_path_coords_nt.csv"):
    coords = pd.read_csv(sys.argv[1] + ".reg_path_coords_nt.csv")
    check_for_benign(sys.argv[1], coords)
else:
    print("No regulated nucleotide regions to clear")
