from process_outputs import *
from utils import *
import sys, os
import pandas as pd

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query name\n")
    exit(1)

biorisk = sys.argv[1] + ".biorisk.hmmsearch"
benign = sys.argv[1] + ".benign.hmmsearch"
synbio = sys.argv[1] + ".benign.blastn"
taxid = sys.argv[1] + ".nr.diamond"

reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

# viz biorisk hits
plot_hmmer(biorisk, 10)

# viz benign hits
plot_hmmer(benign, 10)
plot_blast(file=synbio, nhits = 10, query=sys.argv[1])

# viz taxon IDs
plot_blast(taxid, 10)
plot_tax(taxid, reg_ids, sys.argv[1]) # this is reusing an object already created in taxid screening, so could be accelerated if we generate viz for all queries earlier in the pipeline

