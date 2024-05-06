"""

screen.py

"""

import argparse
import os
import shutil
import subprocess
from common_mechanism.utils import directory_arg, file_arg

def main():
    parser = argparse.ArgumentParser(description='Run pipeline')
    parser.add_argument(dest='fasta_file', type=file_arg, help='FASTA file to screen')
    parser.add_argument(
        '-d',
        '--databases',
        dest='database_dir',
        type=directory_arg,
        required=True,
        help='Path to databases directory'
    )
    parser.add_argument(
        '-o', '--output', dest='output_prefix', help='Output prefix (can be string or directory)'
    )
    parser.add_argument(
        '-t', '--threads', dest='threads', type=int, default=1, help='Threads available'
    )
    parser.add_argument(
        '-p',
        '--protein-search-tool',
        dest='protein_search_tool',
        choices=["blastx", "diamond"],
        default="blastx",
        help='Tool for finding regulated pathogen proteins'
    )
    parser.add_argument(
        '-f',
        '--fast',
        dest='fast_mode',
        action='store_true',
        help='Run in fast mode and skip homology search'
    )
    parser.add_argument(
        '-n',
        '--nt-search',
        dest='nt_search',
        action='store_false',
        help='Skip nucleotide search if no protein hits are found in a region'
    )
    parser.add_argument(
        '-c', '--cleanup', dest='cleanup', action='store_true', help='Delete intermediate files'
    )
    args = parser.parse_args()

    # Get input file
    input_file = args.fasta_file
    input_name = os.path.splitext(os.path.basename(input_file))[0]

    # If output not specified, or only directory, use input filename as  
    if not output_prefix or os.path.isdir(output_prefix):
        output_prefix = output_prefix + input_name

    # Name output file
    screen_file = f"{output_prefix}.screen"

    # Set all the constants
    THREADS = args.threads
    PROTEIN_SEARCH_TOOL = args.protein_search_tool
    IN_FAST_MODE = args.fast_mode
    SKIP_NT_SEARCH = args.nt_search
    DO_CLEANUP = args.cleanup

    # TODO: search input database folder and find the right file locations, rather than requiring hardcoded dir names
    # Check that databases exist

    # Remove spaces or blank characters from the input file
    tmp_output = f"{output_prefix}.tmp"
    cleaned_fasta = f"{output_prefix}.fasta"
    shutil.copyfile(input_file, tmp_output)
    with open(tmp_output, 'r') as fin, open(cleaned_fasta, 'w') as fout:
        for line in fin:
            # Remove whitespace characters (including non-breaking spaces) and replace with underscores
            modified_line = ''.join('_' if c.isspace() or c == '\xc2\xa0' else c for c in line)
            fout.write(modified_line)

    # Biorisk scan
    print(">> STEP 1: Checking for biorisk genes...")
    print("\t...running transeq")

    command = ["transeq", input_file, f"{output_prefix}.transeq.faa", "-frame", "6", "-clean"]
    result = subprocess.run(command, capture_output=True, text=True)

    # Check the return code
    if result.returncode != 0:
        print("\t ERROR: transeq failed")
        with open(screen_file, "a", encoding="utf-8") as f:
            f.write("\t ERROR: transeq failed\n")
    else:
        with open(tmp_output, "a", encoding="utf-8") as f:
            f.write(result.stdout)

if __name__ == "__main__":
    main()
