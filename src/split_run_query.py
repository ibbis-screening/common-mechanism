# function to split multifastas into individual files and screen them

from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser(description="Split the fasta file into individual file with each gene seq")
parser.add_argument('-f', action='store', dest='fasta_file', help='Input fasta file')
result = parser.parse_args()

f_open = open(result.fasta_file, "r")

for rec in SeqIO.parse(f_open, "fasta"):
   id = rec.description
   file_id = "".join(x for x in id if x.isalnum())
   seq = rec.seq
   id_file = open(file_id+".fasta", "w")
   id_file.write(">"+str(id)+"\n"+str(seq))
   id_file.close()
   print(file_id)
   os.system("src/run_pipeline.sh " + file_id + ".fasta > " + file_id + ".screen.txt")

f_open.close()

# usage: python split_query.py -f fasta_file

