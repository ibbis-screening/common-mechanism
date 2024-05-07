"""

screen.py

"""

import argparse
import datetime
import logging
import os
import shutil
import subprocess
from common_mechanism.utils import directory_arg, file_arg

DESCRIPTION = "Run Common Mechanism screening on an input FASTA."

def add_args(parser):
    """
    Add module arguments to an ArgumentParser object.
    """
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
        help='Tool for homology search to identify regulated pathogen proteins'
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
    return parser

def run_as_subprocess(command, out_file, raise_errors=False):
    """
    Run a command using subprocess.run, piping stdout and stderr to `out_file`.
    """
    with open(out_file, "a", encoding='utf-8') as f:
        return subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors)

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()

    # Get input file
    input_file = args.fasta_file
    input_name = os.path.splitext(os.path.basename(input_file))[0]

    # If output not specified, or just directory, use input filename (sans extension) as prefix
    output_prefix = args.output_prefix
    if not output_prefix:
        output_prefix = os.path.splitext(input_file)[0]
    elif os.path.isdir(output_prefix):
        output_prefix = output_prefix + input_name

    # Name output file
    screen_file = f"{output_prefix}.screen"
    tmp_output_file = f"{output_prefix}.tmp"

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(screen_file, 'a')
        ]
    )
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(message)s',
        handlers=[
            logging.FileHandler(tmp_output_file, 'a')
        ]
    )

    # Set all the constants
    THREADS = args.threads
    PROTEIN_SEARCH_TOOL = args.protein_search_tool
    IN_FAST_MODE = args.fast_mode
    SKIP_NT_SEARCH = args.nt_search
    DO_CLEANUP = args.cleanup

    # TODO: search input database folder and find the right file locations,
    # rather than requiring hardcoded dir names
    # Check that databases exist
    DATABASE_DIR = args.database_dir
    BIORISK_DIR = os.path.join(args.database_dir, 'biorisk_db')
    BENIGN_DIR = os.path.join(args.database_dir, 'benign_db')
    if PROTEIN_SEARCH_TOOL == 'blastx':
        NR_DIR = os.path.join(args.database_dir, 'nr_dmnd')
    elif PROTEIN_SEARCH_TOOL == 'diamond':
        NR_DIR = os.path.join(args.database_dir, 'nr_blast')

    # Add a check for transeq as well?
    # check_datbases(BIORISK_DIR, 'biorisk.hmm')
    # check_databases(BENIGN_DIR, 'benign.hmm'

    # Remove spaces or blank characters from the input file
    cleaned_fasta = f"{output_prefix}.fasta"
    shutil.copyfile(input_file, tmp_output_file)
    with (
        open(tmp_output_file, 'r', encoding='utf-8') as fin,
        open(cleaned_fasta, 'w', encoding='utf-8') as fout
    ):
        for line in fin:
            # Remove whitespace (including non-breaking spaces) and replace with underscores
            line = line.strip()
            modified_line = ''.join('_' if c.isspace() or c == '\xc2\xa0' else c for c in line)
            fout.write(f"{modified_line}{os.linesep}")

    # Biorisk scan
    logging.info(">> STEP 1: Checking for biorisk genes...")

    logging.debug("\t...running transeq")
    transeq_out  = f"{output_prefix}.transeq.faa"
    command = ["transeq", input_file, transeq_out, "-frame", "6", "-clean"]
    result = run_as_subprocess(command, tmp_output_file)
    if result.returncode != 0:
        logging.info("\t ERROR: transeq failed")

    logging.debug("\t...running hmmscan")
    biorisk_out = f"{output_prefix}.biorisk.hmmscan"
    command = ["hmmscan", "--domtblout", biorisk_out, f"{BIORISK_DIR}/biorisk.hmm", transeq_out]
    run_as_subprocess(command, tmp_output_file)

    logging.debug("\t...checking hmmscan results")
    command = ["python3", "", biorisk_out, f"{BIORISK_DIR}/biorisk.hmm", transeq_out]
    run_as_subprocess(command, screen_file)

    logging.info(" STEP 1 completed at %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    if IN_FAST_MODE:
        logging.info(" >> FAST MODE: Skipping steps 2-3")
    else:
        logging.info(" >> STEP 2: Checking regulated pathogen proteins...")

if __name__ == "__main__":
    main()
