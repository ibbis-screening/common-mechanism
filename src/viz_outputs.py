import argparse
import re
from os.path import exists
import pandas as pd
from process_outputs import plot_blast, plot_hmmer, plot_blast_frag, plot_tax

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q","--query", dest="query",
        required=True, help="Query name")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="Database folder")
    parser.add_argument("-n","--n-hits", dest="n_hits",
        required=False, help="Number of hits to display")
    args = parser.parse_args()

    biorisk = args.query + ".biorisk.hmmscan"
    benign = args.query + ".benign.hmmscan"
    synbio = args.query + ".benign.blastn"
    taxid = args.query + ".nr.blastx"
    if not exists(taxid):
        taxid = args.query + ".nr.dmnd"
    taxid2 = args.query + ".nt.blastn"

    reg_ids = pd.read_csv(args.db + '/biorisk_db/reg_taxids', header=None)
    vax_ids = pd.read_csv(args.db + '/benign_db/vax_taxids', header=None)
    biorisk_desc = pd.read_csv(args.db + '/biorisk_db/biorisk_annotations.csv')
    benign_desc = pd.read_csv(args.db + '/benign_db/benign_annotations.tsv', sep="\t")

    # define number of hits to display
    if args.n_hits is None:
        args.n_hits = 10
    else:
        args.n_hits = int(args.n_hits)

    # viz biorisk hits
    plot_hmmer(biorisk, biorisk_desc, args.n_hits)

    # viz benign hits
    plot_hmmer(benign, benign_desc, args.n_hits)

    # viz taxon IDs
    plot_blast(taxid, reg_ids=reg_ids, vax_ids=vax_ids, nhits=args.n_hits)
    # if the nc search is on fragments, treat it differently
    if exists(taxid2):
        with open(taxid2) as f:
            txt=f.read()
        if re.findall(':(.+?)-(.+?)$', txt):
            plot_blast_frag(taxid2, reg_ids=reg_ids, vax_ids=vax_ids, nhits=args.n_hits)
        else:
            plot_blast(taxid2, reg_ids=reg_ids, vax_ids=vax_ids, nhits=args.n_hits)

    # this is reusing an object already created in taxid screening, so could be accelerated if we generate viz for all
    # queries earlier in the pipeline
    plot_tax(taxid, reg_ids, args.query)

if __name__ == "__main__":
    main()
