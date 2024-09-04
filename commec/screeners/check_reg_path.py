#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks results for regulated pathogen and prints any matched coordinates. Ignores any
synthetic constructs.

Usage:
  python check_reg_path.py -i INPUT -d database_folder -t threads
      -i, --input 
"""
import logging
import re
import os
import sys
import argparse
import textwrap
import pandas as pd

from commec.utils.file_utils import FileTools
from commec.tools.blast_tools import BlastTools

pd.set_option("display.max_colwidth", 10000)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",dest="in_file", 
        required=True, help="Input query file (e.g. QUERY.nr.dmnd)")
    parser.add_argument("-d","--database", dest="db",
        required=True,help="database folder (must contain vax_taxids and reg_taxids file)")
    parser.add_argument("-t","--threads", dest="threads",
        required=True,help="number of threads")
    args=parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    rv = check_for_regulated_pathogens(args.in_file, args.db, args.threads)
    sys.exit(rv)


def check_for_regulated_pathogens(input_file : str, input_database_dir : str, n_threads : int):
    """ Check an input file (output from a database search) for regulated pathogens, from the benign and biorisk database taxids."""

    #check input files
    if not os.path.exists(input_file):
        logging.error("\t...input query file %s does not exist\n" % input_file)
        return 1
    if not os.path.exists(input_database_dir + "/benign_db/vax_taxids.txt"):
        logging.error("\t...benign db file %s does not exist\n" % (input_database_dir + "/benign_db/vax_taxids.txt"))
        return 1
    if not os.path.exists(input_database_dir + "/biorisk_db/reg_taxids.txt"):
        logging.error("\t...biorisk db file %s does not exist\n" % (input_database_dir + "/biorisk_db/reg_taxids.txt"))
        return 1

    # read in files
    reg_ids = pd.read_csv(input_database_dir + "/biorisk_db/reg_taxids.txt", header=None)
    vax_ids = pd.read_csv(input_database_dir + "/benign_db/vax_taxids.txt", header=None)

    sample_name = re.sub(".nr.*", "", input_file)
    sample_name = re.sub(".nt.blastn", "", sample_name)
    
    # if there are already regulated regions written to file for this query, add to them
    hits1 = None
    if os.path.exists(sample_name + ".reg_path_coords.csv"):
        hits1 = pd.read_csv(sample_name + ".reg_path_coords.csv")

    if FileTools.is_empty(input_file):
        logging.info("\tERROR: Homology search has failed\n")
        return 1

    if not FileTools.has_hits(input_file):
        logging.info("\t...no hits\n")
        return 0

    blast = BlastTools.readblast(input_file)
    blast = BlastTools.taxdist(blast, reg_ids, vax_ids, input_database_dir + "/taxonomy/", n_threads)
    blast = blast[blast['species'] != ""] # ignore submissions made above the species level

    # trim down to the top hit for each region, ignoring any top hits that are synthetic constructs
    #interesting_cols = ['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']

    blast2 = BlastTools.trimblast(blast)
    blast2 = BlastTools.tophits(blast2) # trims down to only label each base with the top matching hit, but includes the different taxids attributed to the same hit

    reg_bac = 0
    reg_vir = 0
    reg_fung = 0

    # if this is the nucleotide screen, check if any weak protein flags can be negated with strong non-regulated nt ones
    if re.findall(".nt.blastn", input_file):
        if hits1 is not None:
            for region in range(0, hits1.shape[0]): # for each regulated pathogen region
                # look at only the hits that overlap it
                htrim = blast2[~((blast2['q. start'] > hits1['q. end'][region]) & (blast2['q. end'] > hits1['q. end'][region])) & ~((blast2['q. start'] < hits1['q. start'][region]) & (blast2['q. end'] < hits1['q. start'][region]))]
                species_list = textwrap.fill(", ".join(set(htrim['species'])), 100).replace("\n", "\n\t\t     ")
                taxid_list = textwrap.fill(", ".join(map(str, set(htrim['subject tax ids']))), 100).replace("\n", "\n\t\t     ")
                percent_ids = (" ".join(map(str, set(htrim['% identity']))))
                if htrim.shape[0] > 0:
                    if any(htrim['q. coverage'] > 0.90):
                        htrim = htrim[htrim['q. coverage'] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim['subject title'][row]
                            descriptions.append(hit)
                        annot_string = "\n\t...".join(str(v) for v in descriptions)
                        logging.info("\t...Regulated protein region at bases " + str(int(hits1['q. start'][region])) + " to " + str(int(hits1['q. end'][region])) + " overlapped with a nucleotide hit\n")
                        logging.info("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))


    if blast2['regulated'].sum(): # if ANY of the trimmed hits are regulated
        hits = pd.DataFrame(columns=['q. start', 'q. end'])
        with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        # for each hit (subject acc) linked with at least one regulated taxid
            for site in set(blast2['q. start'][blast2['regulated'] != False]):
                subset = blast2[(blast2['q. start'] == site)]
                subset = subset.sort_values(by=['regulated'], ascending=False)
                subset = subset.reset_index(drop=True)
                org = ""
            
                blast2 = blast2.dropna(subset = ['species'])
                n_reg = (blast2['regulated'][blast2['q. start'] == site] != False).sum()
                n_total = len(blast2['regulated'][blast2['q. start'] == site])
                gene_names = ", ".join(set(subset['subject acc.']))
                end = blast2['q. end'][blast2['q. start'] == site].max()
                coordinates = str(int(site)) + " - " + str(int(end))

                species_list = textwrap.fill(", ".join(set(blast2['species'][blast2['q. start'] == site])), 100).replace("\n", "\n\t\t     ")
                desc = blast2['subject title'][blast2['q. start'] == site].values[0]
                taxid_list = textwrap.fill(", ".join(map(str, set(blast2['subject tax ids'][blast2['q. start'] == site]))), 100).replace("\n", "\n\t\t     ")
                percent_ids = (" ".join(map(str, set(blast2['% identity'][blast2['q. start'] == site]))))
                reg_ids = (" ".join(map(str, set(blast2['regulated'][(blast2['q. start'] == site) & (blast2['regulated'] != False)]))))

                # if some of the organisms with this sequence aren't regulated, say so
                if (n_reg < n_total):
                    logging.info("\t\t --> Best match to sequence(s) %s at bases %s found in both regulated and non-regulated organisms\n" % (gene_names, coordinates))
                    logging.info("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))
                    logging.info("\t\t     Description: %s\n" % (desc))
                    # could explicitly list which are and aren't regulated?
                # otherwise, raise a flag and say which superkingdom the flag belongs to
                elif (n_reg == n_total):
                    if subset['superkingdom'][0] == "Viruses":
                        reg_vir = 1
                        org = "virus"
                    elif subset['superkingdom'][0] == "Bacteria": 
                        reg_bac = 1
                        org = "bacteria"
                    elif 'superkingdom' in subset:
                        if subset['superkingdom'][0] == "Eukaryota":
                            org = "eukaryote"
                            reg_fung = 1

                    new_hits = subset[['q. start', 'q. end']].dropna()
                    if not new_hits.empty and not hits.empty:
                         hits = pd.concat([hits, new_hits], ignore_index=True)
                    elif not new_hits.empty:
                        hits = new_hits.copy()
                    logging.info("\t\t --> Best match to sequence(s) %s at bases %s found in only regulated organisms: FLAG (%s)\n" % (gene_names, coordinates, org))
                    logging.info("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))
                    logging.info("\t\t     Description: %s\n" % (desc))
                else: # something is wrong, n_reg > n_total
                    logging.info("\t...gene: %s\n" % gene_names)
                    logging.info("%s\n" % (blast['regulated'][blast['subject acc.'] == gene_names]))
        hits = hits.drop_duplicates()
        #Create output file 
        if hits1 is not None:
            hits = pd.concat([hits1, hits])
        hits.to_csv(sample_name + ".reg_path_coords.csv", index=False)

    if reg_vir == 0 and reg_bac == 0 and reg_fung == 0 and reg_fung == 0:
        logging.info("\t\t --> no top hit exclusive to a regulated pathogen: PASS\n")

    return 0

if __name__ == "__main__":
    main()
