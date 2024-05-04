"""
function to split multifastas into individual files and screen them
usage: python split_query.py -f fasta_file
"""

import argparse
import string
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Split the fasta file into individual file with each gene seq")
parser.add_argument('-f', action='store', dest='fasta_file', help='Input fasta file')
result = parser.parse_args()

with open(result.fasta_file, "r", encoding="utf-8") as input_file:
    for record in SeqIO.parse(input_file, "fasta"):
        desc = record.description
        desc = desc.strip()

        VALID_FILENAME_CHARS = f"-.{string.ascii_letters}{string.digits}"
        FILENAME = "".join(x for x in desc if x in VALID_FILENAME_CHARS)
        if len(FILENAME) > 150:
            FILENAME = FILENAME[:150]

        with open(f"{FILENAME}.fasta", "w", encoding="utf-8") as output_file:
            output_file.write(">"+str(desc)+"\n"+str(record.seq))
