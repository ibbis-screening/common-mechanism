"""
Split a FASTA file containing multiple records into individual files, one for each record.

You can call it as a script:
    
    python split.py input.fasta 
"""
import argparse
import os
import string
from Bio import SeqIO

VALID_FILENAME_CHARS = f"-._{string.ascii_letters}{string.digits}"

def clean_description(description):
    """
    Cleans the description from a sequence record for use as part of a filename.
    """
    cleaned = description.strip()
    cleaned = "".join(x for x in cleaned if x in VALID_FILENAME_CHARS)
    if len(cleaned) > 150:
        cleaned = cleaned[:150]
    return cleaned

def write_split_fasta(fasta_file):
    """
    Parse all sequence records in an input FASTA file, and write a new file for each record.
    """
    output_dir = os.path.dirname(fasta_file)

    with open(fasta_file, "r", encoding="utf-8") as input_file:
        for record in SeqIO.parse(input_file, "fasta"):
            desc = clean_description(record.description)
            outpath = os.path.join(output_dir, f"{desc}.fasta")

            # Avoid overwriting input files
            if outpath == fasta_file:
                outpath = os.path.join(output_dir, f"{desc}-split.fasta")

            with open(outpath, "w", encoding="utf-8") as output_file:
                output_file.write(f">{desc}{os.linesep}")
                output_file.write(record.seq)

def main():
    """
    Main function. Passes FASTA file to `write_split_fasta`.

    Arguments:
      - fasta_file: Path to the input FASTA file.

    """
    parser = argparse.ArgumentParser(
        description="Split a FASTA file containing multiple records into individual files, one for each record."
    )
    parser.add_argument(action='store', dest='fasta_file', help='Input fasta file')
    result = parser.parse_args()

    write_split_fasta(result.fasta_file)

if __name__ == "__main__":
    main()
