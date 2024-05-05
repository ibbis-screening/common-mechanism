"""
Split FASTA files containing multiple records into individual files, one for each record.

usage: python split.py input.fasta 
"""

import argparse
import os
import string
from Bio import SeqIO

VALID_FILENAME_CHARS = f"-.{string.ascii_letters}{string.digits}"

def clean_description(description):
    description = description.strip()
    record_name = "".join(x for x in description if x in VALID_FILENAME_CHARS)
    if len(record_name) > 150:
        record_name = record_name[:150]

    return f"{record_name}.fasta"


def write_split_fasta(fasta_file):
    output_dir = os.path.dirname(fasta_file)
    with open(fasta_file, "r", encoding="utf-8") as input_file:
        for record in SeqIO.parse(input_file, "fasta"):
            desc = clean_description(record.description)
            out_filename = os.path.join(output_dir, f"{desc}.fasta")
            with open(out_filename, "w", encoding="utf-8") as output_file:
                print(f">{desc}", file=output_file)
                print(record.seq, file=output_file)

def main():
    parser = argparse.ArgumentParser(
        description="Split a FASTA file containing multiple records into individual files, one for each record."
    )
    parser.add_argument(action='store', dest='fasta_file', help='Input fasta file')
    result = parser.parse_args()

    write_split_fasta(result.fasta_file)

if __name__ == "__main__":
    main()
