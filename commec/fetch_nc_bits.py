#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Fetch parts of a query that had no high-quality protein matches for use in nucloetide screening.

Usage:
    fetch_nc_bits.py query_name fasta_file_path
"""
import sys
import shutil
import re
from Bio import SeqIO
from commec.utils import readblast, trimblast, is_empty, has_hits


def get_high_identity_matches(blast_output_file, threshold=90):
    """Read all hits with high sequence identity from a BLAST results file."""
    hits = readblast(blast_output_file)
    hits = trimblast(hits)
    return hits[hits["% identity"] >= threshold]


def get_ranges_with_no_hits(blast_df):
    """
    Get indices not covered by the query start / end ranges in the BLAST results.
    """
    hits = blast_df[["q. start", "q. end"]].values.tolist()
    print(hits)
    # Sort each pair to ensure that start < end, then sort entire list of ranges
    hits = sorted([sorted(pair) for pair in hits])

    nc_ranges = []

    # Include the start if the first hit begins more than 50 bp after the start
    if hits[0][0] > 50:
        nc_ranges.append([1, hits[0][0] - 1])

    # Add ranges if there is a noncoding region of >=50 between hits
    for i in range(len(hits) - 1):
        nc_start = hits[i][1] + 1      # starts after this hit
        nc_end = hits[i + 1][0] - 1    # ends before next hit

        if nc_end - nc_start + 1 >= 50:
            nc_ranges.append([nc_start, nc_end])

    # Include the end if the last hit ends more than 50 bp before the end
    query_length = blast_df['query length'][0]
    if query_length - hits[-1][1] >= 50:
        nc_ranges.append([hits[-1][1] + 1, int(query_length)])

    return nc_ranges


def get_records(fasta_file_path):
    """Parse SeqIO records from input FASTA."""
    with open(fasta_file_path, "r", encoding="utf-8") as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        return records


def write_nc_sequences(nc_ranges, record, outfile: str):
    """
    Write a FASTA containing only the parts of the record without a high-quality protein match.
    """
    nc_sequences = []
    for start, stop in nc_ranges:
        # subtract 1 just from `start` to adjust to 0-based index and capture full range
        sequence = record.seq[int(start) - 1 : int(stop)]
        nc_sequences.append(f">{record.id} {record.description} {start}-{stop}\n{sequence}\n")

    with open(outfile, "w", encoding="utf-8") as output_file:
        output_file.writelines(nc_sequences)


def main(protein_results, query_fasta):
    """Fetch noncoding regions > 50bp and write to a new file."""
    outfile = re.sub(".nr.*", "", protein_results) + ".noncoding.fasta"

    if is_empty(protein_results) or not has_hits(protein_results):
        sys.stdout.write("\t...no protein hits found, screening entire sequence\n")
        shutil.copyfile(query_fasta, outfile)
        return

    blast_df = get_high_identity_matches(protein_results)

    if blast_df.empty:
        sys.stdout.write(
            "\t...protein hits all low percent identity (<90%) - screening entire sequence\n"
        )
        shutil.copyfile(query_fasta, outfile)
        return

    sys.stdout.write(
        "\t...protein hits found, fetching nt regions not covered by a 90% ID hit or better\n"
    )
    ranges_to_screen = get_ranges_with_no_hits(blast_df)

    # if the entire sequence, save regions <50 bases, is covered with protein, skip nt scan
    if not ranges_to_screen:
        sys.stdout.write("\t\t --> no noncoding regions >= 50 bases found, skipping nt scan\n")
        return

    records = get_records(query_fasta)

    if len(records) > 1:
        sys.stdout.write(
            "\t\t...WARNING: Only fetching nucleotides from first record in multifasta: "
            + f"{query_fasta}\n"
        )

    sys.stdout.write(
        "\t\t Fetching the following noncoding regions: "
        + ", ".join(f"{start}-{end}" for start, end in ranges_to_screen)
        + "\n"
    )
    write_nc_sequences(ranges_to_screen, records[0], outfile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    protein_results, query_fasta = sys.argv[1], sys.argv[2]
    main(protein_results, query_fasta)
