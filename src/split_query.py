# function to split multifastas into individual files and screen them
# usage: python split_query.py -f fasta_file

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Split the fasta file into individual file with each gene seq")
parser.add_argument('-f', action='store', dest='fasta_file', help='Input fasta file')
result = parser.parse_args()

f_open = open(result.fasta_file, "r")

for rec in SeqIO.parse(f_open, "fasta"):
    id = rec.description
    file_id = "".join(x for x in id if x.isalnum())
    if len(file_id) > 150:
        file_id = file_id[:150]
    seq = rec.seq
    id_file = open(file_id+".fasta", "w")
    id_file.write(">"+str(id)+"\n"+str(seq))
    id_file.close()
    print(file_id)
f_open.close()
