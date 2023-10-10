# usage: fetch_nc_bits.py query_name fasta_file_path

# checks whether there are any hits to nr for a query. If there aren't any over a gien significance level, prints the whole sequence to a noncoding query file. If there are hits, fetches the nucleotide regions between these hits and singles them out for nucleotide screening

from utils import *
import sys, shutil
# import pybedtools
import re
from Bio import SeqIO

query = sys.argv[1]
f_file = sys.argv[2]

nc_bits = 0

# check if the nr hits file is empty
if check_blastfile(query) == 0:
    nc_bits = "all"
elif check_blastfile(query) == 2:
    sys.stdout.write("\t...no hits to the nr database\n")
    nc_bits = "all"
# if not, check whether any of the hits has an E-value > 1e-30
# if so, find the start and end of those hits and use these to get the coordinates of non-coding regions
else:
    blast = readblast(query)
    blast = trimblast(blast)
    # print(blast['% identity'])
    blast = blast[blast['% identity'] >= 90]
    if blast.shape[0] > 0:
    # find noncoding bits
        sys.stdout.write("\t...protein hits found, fetching nt regions not covered by a 90% ID hit or better\n")
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
        sys.stdout.write("\t...protein hits all low percent identity (<90%) - screening entire sequence\n")
        nc_bits = "all"

# fetch noncoding sequences

outfile = re.sub(".nr.*", "", query) + '.noncoding.fasta'
# print("Outfile: " + outfile)

def fetch_sequences(seqid, nc_bits, f_file, outfile):
    tofetch = []
    for (start, stop) in nc_bits:
        tofetch.append((seqid, start, stop))

    if tofetch:
        with open(f_file, "r") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            sequences = []
            print(tofetch)
            for (seqid, start, stop) in tofetch:
                start = int(start)
                stop = int(stop)
                for record in records:
                    print("$" + str(seqid) + "$ $" + str(record.id))
                    if str(record.id) == str(seqid):
                        sequence = record.seq[start - 1 : stop]  # Adjust start to 0-based index
                        print(sequence, record.seq)
                        sequences.append(f">{seqid} {start}-{stop}\n{sequence}\n")
                        break
                    else:
                        print("Error")
        print(sequences)
        with open(outfile, "w") as output_file:
            output_file.writelines(sequences)

if nc_bits == "all":
    shutil.copyfile(f_file, outfile)
    # print("\t...no significant protein hits")
elif nc_bits == []: # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
    sys.stdout.write("\t\t --> no noncoding regions >= 50 bases found, skipping nt scan\n")
else: 
    print("pulling out noncoding bits")
    seqid = blast.iloc[0][0]
    fetch_sequences(seqid, nc_bits, f_file, outfile)






