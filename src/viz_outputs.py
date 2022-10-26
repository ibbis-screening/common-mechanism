from process_outputs import *
from utils import *
import sys, os
import pandas as pd
from os.path import exists

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query name\n")
    exit(1)

biorisk = sys.argv[1] + ".biorisk.hmmsearch"
benign = sys.argv[1] + ".benign.hmmsearch"
synbio = sys.argv[1] + ".benign.blastn"
taxid = sys.argv[1] + ".nr.blastx"
if not exists(taxid):
    taxid = sys.argv[1] + ".nr.dmnd"
taxid2 = sys.argv[1] + ".nt.blastn"

reg_ids = pd.read_csv(os.environ['DB_PATH'] + '/biorisk/reg_taxids', header=None)

# viz biorisk hits
plot_hmmer(biorisk, 10)

# viz benign hits
plot_hmmer(benign, 10)
plot_blast(file=synbio, nhits = 10, query=sys.argv[1])

# viz taxon IDs
plot_blast(taxid, 10)
plot_blast(taxid2, 10)
plot_tax(taxid, reg_ids, sys.argv[1]) # this is reusing an object already created in taxid screening, so could be accelerated if we generate viz for all queries earlier in the pipeline

