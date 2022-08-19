# usage: fetch_nc_bits.py query_name fasta_file_path

from utils import *
import os, sys
import pandas as pd
import taxoniq
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from collections import defaultdict
import pybedtools

#Input parameter error checking
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)
    
query = sys.argv[1]
f_file = sys.argv[2]
blast_file = query + ".nr.blastx"

# find noncoding bits
blast = readblast(blast_file)
blast = trimblast(blast)
# print(diamond)
hits = []
for i in range(blast.shape[0]):
    pair = [blast['q. start'][i], blast['q. end'][i]]
    pair.sort
    # print(pair)
    hits.append(pair)
#print(hits)

nc_bits = []
if hits[0][0] != 1:
    nc_bits.append([1,hits[0][0]])
for i in range(len(hits)-1):
    if hits[i][1] < hits[i+1][0]:
        nc_bits.append([hits[i][1], hits[i+1][0]])
#print(nc_bits)


# fetch noncoding sequences

query="1_biorisk"
seqid = blast.iloc[0][0]

# nc_bits.append([2,4])

tofetch = ""
for (start, stop) in nc_bits:
    tofetch = tofetch + str(seqid) + " " + str(start) + " " + str(stop) + "\n"

a = pybedtools.BedTool(tofetch, from_string=True)
fasta = f_file
outfile = query + '_nc.fasta'
a = a.sequence(fi=fasta, fo=outfile)
# print(open(a.seqfn).read())


