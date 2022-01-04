from utils import *
import os, sys
import pandas as pd
#import SeqIO

#Input parameter error checking 
if len(sys.argv) < 2:
    sys.stdout.write("\tERROR: Please provide two files: one sequence name, one FASTA for the sequence\n")
    exit(1)

if os.path.exists(sys.argv[1] + ".reg_path_coords.csv"):
    coords = pd.read_csv(sys.argv[1] + ".reg_path_coords.csv")
    
    # for each set of hits, need to pull out the coordinates covered by benign entries
    # overall the majority (?) of regulated pathogen coordinates must be covered with benign content
    hmmsearch = sys.argv[1] + ".benign.hmmsearch"
    qlen = 0
    if checkfile(hmmsearch) == 2:
        print("No housekeeping genes found")
        exit(1)
    
    hmmer = readhmmer(hmmsearch)
    qlen = hmmer['qlen'][0]
    for region in range(0, coords.shape[0]): # for each regulated pathogen region
        # look at only the hmmer hits that overlap with it
        htrim = hmmer[~((hmmer['ali from'] > coords['q. end'][region]) & (hmmer['ali to'] > coords['q. end'][region])) & ~((hmmer['ali from'] < coords['q. start'][region]) & (hmmer['ali to'] < coords['q. start'][region]))]
        if htrim.shape[0] > 0:
            htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'][0])
#            htrim['coverage'] = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'][0]
            if any(htrim['coverage'] > 0.90):
                print("Housekeeping genes - >90% coverage achieved")
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
            htrim = blastn[~((blastn['q.start'] > coords['q. end'][region]) & (blastn['q. end'] > coords['q. end'][region])) & ~((blastn['q.start'] < coords['q. start'][region]) & (blastn['q. from'] < coords['q. start'][region]))]
            htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'][0])
            if any(htrim2['coverage'] > 0.90):
                print("Housekeeping genes - >90% coverage achieved - PASS")
else:
    print("No regulated regions to clear")
