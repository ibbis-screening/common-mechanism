# usage: fetch_nc_bits.py query_name fasta_file_path

# checks whether there are any hits to nr for a query. If there aren't any over a gien significance level, prints the whole sequence to a noncoding query file. If there are hits, fetches the nucleotide regions between these hits and singles them out for nucleotide screening

from utils import *
import sys, shutil
import pybedtools

#Input parameter error checking
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)
    
query = sys.argv[1]
f_file = sys.argv[2]
outname = query.replace(".nr.*", "")

# check if the nr hits file is empty
if check_blastfile(query) == 0 or check_blastfile(query) == 2:
    nc_bits = "all"
# if not, check whether any of the hits has an E-value > 1e-30
# if so, find the start and end of those hits and use these to get the coordinates of non-coding regions
else:
    blast = readblast(query)
    blast = trimblast(blast)
    blast = blast[blast['evalue'] < 1e-30]
    if blast.shape[0] > 0:
    # find noncoding bits
        hits = []
        for i in range(blast.shape[0]):
            pair = [blast['q. start'][i], blast['q. end'][i]]
            pair.sort()
            hits.append(pair)
        hits = sorted(hits, key=lambda x: x[0])

        nc_bits = []
        print(hits)
        if hits[0][0] >50:
            nc_bits.append([1,hits[0][0]])
        for i in range(len(hits)-1):
            if hits[i][1] < (hits[i+1][0] - 49): # if there's a noncoding region of >=50 between hits
                nc_bits.append([hits[i][1], hits[i+1][0]])
    else:
        nc_bits = "all"

# fetch noncoding sequences

outfile = query + '_nc.fasta'
if nc_bits == "all":
    shutil.copyfile(f_file, outfile)
    print("no significant protein hits")
else: 
    print("pulling out noncoding bits")
    seqid = blast.iloc[0][0]

    print(nc_bits)

    tofetch = ""
    for (start, stop) in nc_bits:
        tofetch = tofetch + str(seqid) + " " + str(start) + " " + str(stop) + "\n"

    if tofetch != "":
        a = pybedtools.BedTool(tofetch, from_string=True)
        print(a)
        fasta = f_file
        a = a.sequence(fi=fasta, fo=outfile)


