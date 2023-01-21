# usage: fetch_nc_bits.py query_name fasta_file_path

# checks whether there are any hits to nr for a query. If there aren't any over a gien significance level, prints the whole sequence to a noncoding query file. If there are hits, fetches the nucleotide regions between these hits and singles them out for nucleotide screening

from utils import *
import sys, shutil
import pybedtools
import re

query = sys.argv[1]
f_file = sys.argv[2]

nc_bits = 0

# check if the nr hits file is empty
if check_blastfile(query) == 0:
    sys.stdout.write("\tERROR: Protein search has failed\n")
    nc_bits = "all"
elif check_blastfile(query) == 2:
    sys.stdout.write("\t...no hits to the nr database\n")
    nc_bits = "all"
# if not, check whether any of the hits has an E-value > 1e-30
# if so, find the start and end of those hits and use these to get the coordinates of non-coding regions
else:
    blast = readblast(query)
    blast = trimblast(blast)
    blast = blast[blast['evalue'] < 1e-30]
    # print(blast)
    if blast.shape[0] > 0:
    # find noncoding bits
        sys.stdout.write("\t...protein hits found, fetching noncoding regions\n")
        hits = []
        for i in range(blast.shape[0]):
            pair = [blast['q. start'][i], blast['q. end'][i]]
            pair.sort()
            hits.append(pair)
        hits = sorted(hits, key=lambda x: x[0])

        # print(hits)
        nc_bits = []
        if hits[0][0] >50:
            nc_bits.append([1,hits[0][0]])
        for i in range(len(hits)-1):
            if hits[i][1] < (hits[i+1][0] - 49): # if there's a noncoding region of >=50 between hits
                nc_bits.append([hits[i][1], hits[i+1][0]])
    else:
        sys.stdout.write("\t...protein hits all low significance - screening entire sequence\n")
        nc_bits = "all"

# fetch noncoding sequences

# print("Outfile: " + outfile)
outfile = re.sub(".nr.*", "", query) + '.noncoding.fasta'
# print(outfile)

if nc_bits == "all":
    shutil.copyfile(f_file, outfile)
    print("\t...no significant protein hits")
elif nc_bits == []: # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
    sys.stdout.write("\t...no noncoding regions >= 50 bases found, skipping nt scan\n")
else: 
    print("pulling out noncoding bits")
    seqid = blast.iloc[0][0]

    tofetch = ""
    for (start, stop) in nc_bits:
        tofetch = tofetch + str(seqid) + " " + str(start) + " " + str(stop) + "\n"

    if tofetch != "":
        a = pybedtools.BedTool(tofetch, from_string=True)
        fasta = f_file
        a = a.sequence(fi=fasta, fo=outfile)
    else:
        print("Error")
        print(nc_bits)


