#! /usr/bin/env python

##################################################################################################
#check_reg_path.py checks results for regulated pathogen and prints any 
#matched coordinates
#   This file will ignore any synthetic constructs
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
##################################################################################################
#Usage:
#   python check_reg_path.py -i INPUT -d database_folder -t threads
#       -i, --input 
##################################################################################################
from utils import *
import os, sys, argparse
import pandas as pd
import textwrap

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
    
    #check input files
    if (not os.path.exists(args.in_file)):
        sys.stderr.write("\t...input query file %s does not exist\n" % args.in_file)
        exit(1)
    if (not os.path.exists(args.db + "/benign_db/vax_taxids")):
        sys.stderr.write("\t...benign db file %s does not exist\n" % (args.db + "/benign_db/vax_taxids"))
        exit(1)
    if (not os.path.exists(args.db + "/biorisk_db/reg_taxids")):
        sys.stderr.write("\t...biorisk db file %s does not exist\n" % (args.db + "/biorisk_db/reg_taxids"))
        exit(1)
    # read in files
    reg_ids = pd.read_csv(args.db + "/biorisk_db/reg_taxids", header=None)
    vax_ids = pd.read_csv(args.db + "/benign_db/vax_taxids", header=None)

    # sample_name = re.sub("\..*", "", args.in_file)
    sample_name = re.sub(".nr.*", "", args.in_file)
    sample_name = re.sub(".nt.blastn", "", sample_name)
    # print("Sample name: " + sample_name)
    
    # if there are already regulated regions written to file for this query, add to them
    hits1 = None
    if os.path.exists(sample_name + ".reg_path_coords.csv"):
        hits1 = pd.read_csv(sample_name + ".reg_path_coords.csv")

    if check_blastfile(args.in_file) == 0:
        sys.stdout.write("\tERROR: Homology search has failed\n")
        exit(1)
    if check_blastfile(args.in_file) == 2:
        sys.stdout.write("\t...no hits\n")
        exit(1)
    blast = readblast(args.in_file)                  #function in utils.py
    blast = taxdist(blast, reg_ids, vax_ids, args.db, args.threads) #function in utils.py
    # print(blast['subject tax ids'])
    blast = blast[(blast['superkingdom'] != "Bacteria") | (blast['species'] != "")] # ignore submissions made above the species level

    # trim down to the top hit for each region, ignoring any top hits that are synthetic constructs
    # print("Original blast")
    # print(blast[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']].sort_values('q. start').head(40))
    blast2 = trimblast(blast)
    # print("Trim blast")
    # print(blast2[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']].sort_values('q. start').head(40))
    blast2 = tophits(blast2) # trims down to only label each base with the top matching hit, but includes the different taxids attributed to the same hit
    # print("Top blast")
    # print(blast2[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']].sort_values('q. start').head(40))

    reg_bac = 0
    reg_vir = 0
    reg_fung = 0

    # if this is the nucleotide screen, check if any weak protein flags can be negated with strong non-regulated nt ones
    if re.findall(".nt.blastn", args.in_file):
        if hits1 is not None:
            for region in range(0, hits1.shape[0]): # for each regulated pathogen region
                # look at only the hits that overlap it
                htrim = blast2[~((blast2['q. start'] > hits1['q. end'][region]) & (blast2['q. end'] > hits1['q. end'][region])) & ~((blast2['q. start'] < hits1['q. start'][region]) & (blast2['q. end'] < hits1['q. start'][region]))]
                species_list = textwrap.fill(", ".join(set(htrim['species'])), 100).replace("\n", "\n\t\t     ")
                taxid_list = textwrap.fill(", ".join(map(str, set(htrim['subject tax ids']))), 100).replace("\n", "\n\t\t     ")
                percent_ids = (" ".join(map(str, set(htrim['% identity']))))
                # print(htrim)
                if htrim.shape[0] > 0:
                    if any(htrim['q. coverage'] > 0.90):
                        htrim = htrim[htrim['q. coverage'] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim['subject title'][row]
                            descriptions.append(hit)
                        annot_string = "\n\t...".join(str(v) for v in descriptions)
                        sys.stdout.write("\t...Regulated protein region at bases " + str(int(hits1['q. start'][region])) + " to " + str(int(hits1['q. end'][region])) + " overlapped with a nucleotide hit\n")
                        sys.stdout.write("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))


    if blast2['regulated'].sum(): # if ANY of the trimmed hits are regulated
        hits = pd.DataFrame(columns=['q. start', 'q. end'])
        # print(blast2)
        # sys.stdout.write("\t...regulated pathogen sequence: PRESENT\n")
        with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
            #  print(blast[['query acc.', 'subject acc.', 'regulated', 'species', 'q. start', 'q. end', '% identity']])
        # for each hit (subject acc) linked with at least one regulated taxid
            for site in set(blast2['q. start'][blast2['regulated'] != False]):
                subset = blast2[(blast2['q. start'] == site)]
                subset = subset.sort_values(by=['regulated'], ascending=False)
                subset = subset.reset_index(drop=True)
                org = ""
            
                # if sum(subset['regulated']) > 0: # if there's at least one regulated hit
                # print("screening at " + str(site))
                # print(subset[['subject acc.', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']])
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
                    sys.stdout.write("\t\t --> Best match to sequence(s) %s at bases %s found in both regulated and non-regulated organisms\n" % (gene_names, coordinates))
                    sys.stdout.write("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))
                    sys.stdout.write("\t\t     Description: %s\n" % (desc))
                    # could explicitly list which are and aren't regulated?
                    # otherwise, raise a flag and say which superkingdom the flag belongs to
                elif (n_reg == n_total):
                    if subset['superkingdom'][0] == "Viruses":
                        reg_vir = 1
                        org = "virus"
                    elif subset['superkingdom'][0] == "Bacteria": 
                        reg_bac = 1
                        org = "bacteria"
                    elif 'kingdom' in subset:
                        if subset['kingdom'][0] == "Fungi":
                            org = "fungi"
                            reg_fung = 1
                    if subset['phylum'][0] == "Oomycota":
                        org = "oomycete"
                        reg_fung = 1 # sorry! to save complexity
                    # sys.stdout.write("\t...%s\n" % (subset['superkingdom'][0]))
                    hits = pd.concat([hits, subset[['q. start', 'q. end']]])
                    sys.stdout.write("\t\t --> Best match to sequence(s) %s at bases %s found in only regulated organisms: FLAG (%s)\n" % (gene_names, coordinates, org))
                    sys.stdout.write("\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n" % (species_list, taxid_list, percent_ids))
                    sys.stdout.write("\t\t     Regulated taxid(s): %s\n" % (reg_ids))
                    sys.stdout.write("\t\t     Description: %s\n" % (desc))
                else: # something is wrong, n_reg > n_total
                    sys.stdout.write("\t...gene: %s\n" % gene_names)
                    sys.stdout.write("%s\n" % (blast['regulated'][blast['subject acc.'] == gene_names]))
        # hits = blast2[blast2['regulated']==True][['q. start', 'q. end']]  # print out the start and end coordinates of the query sequence
        hits = hits.drop_duplicates()
        #Create output file 
        if hits1 is not None:
            hits = pd.concat([hits1, hits])
        # print(hits)
        hits.to_csv(sample_name + ".reg_path_coords.csv", index=False)

    if reg_vir == 0 and reg_bac == 0 and reg_fung == 0 and reg_fung == 0:
        sys.stdout.write("\t\t --> no top hit exclusive to a regulated pathogen: PASS\n")

if __name__ == "__main__":
    main()
