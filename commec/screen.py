#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Run Common Mechanism screening on an input FASTA.

screen.py
        -d DATBASE_FOLDER
        INPUT_FASTA
    Optional parameters:
        -t THREADS (default: 1)
        -o OUTPUT (output prefix, default: input name)
        -p SEARCH_TOOL (default: blastx) 
        -f = use fast mode (default: false)
        -n = skip nucleotide search if no protein hits are found in a region (default: false)
        -c = clean up intermediate files (default: false)

Command-line usage:
    screen.py -d /path/to/databases input.fasta -t 1 -o out_prefix 

"""
import argparse
import datetime
import glob
import logging
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from commec.utils import directory_arg, file_arg, is_likely_protein 
from commec.screen_databases import ScreenDatabases

DESCRIPTION = "Run Common Mechanism screening on an input FASTA."

def add_args(parser):
    """
    Add module arguments to an ArgumentParser object.
    """
    parser.add_argument(dest="fasta_file", type=file_arg, help="FASTA file to screen")
    parser.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        required=True,
        help="Path to databases directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_prefix",
        help="Output prefix (can be string or directory)",
    )
    parser.add_argument(
        "-t", "--threads", dest="threads", type=int, default=1, help="Threads available"
    )
    parser.add_argument(
        "-p",
        "--protein-search-tool",
        dest="protein_search_tool",
        choices=["blastx", "diamond"],
        default="blastx",
        help="Tool for homology search to identify regulated pathogen proteins",
    )
    parser.add_argument(
        "-f",
        "--fast",
        dest="fast_mode",
        action="store_true",
        help="Run in fast mode and skip homology search",
    )
    parser.add_argument(
        "-n",
        "--skip-nt",
        dest="skip_nt_search",
        action="store_true",
        help="Skip nucleotide search if no protein hits are found in a region",
    )
    parser.add_argument(
        "-c",
        "--cleanup",
        dest="cleanup",
        action="store_true",
        help="Delete intermediate files",
    )
    return parser


@dataclass
class ScreenParameters:
    """
    Parameters used for screening; provided by parsing arguments. All have default values.
    """
    threads: int
    protein_search_tool: str
    in_fast_mode: bool
    skip_nt_search: bool
    do_cleanup: bool

def get_output_prefix(input_file, prefix_arg=""):
    """
    Returns a prefix that can be used for all output files.
    
    - If no prefix was given, use the input filename.
    - If a directory was given, use the input filename as file prefix within that directory.
    """
    if not prefix_arg:
        return os.path.splitext(input_file)[0]
    if prefix_arg.endswith("/") or prefix_arg in {".", ".."} or prefix_arg.endswith("\\"):
        # Make the directory if it doesn't exist
        if not os.path.isdir(prefix_arg):
            os.makedirs(prefix_arg)
        # Use the input filename as a prefix within that directory (stripping out the path)
        input_name = os.path.splitext(os.path.basename(input_file))[0]
        return prefix_arg + input_name
    # Existing, non-path prefixes can be used as-is
    return prefix_arg


def get_cleaned_fasta(input_file, out_prefix):
    """
    Return a FASTA where whitespace (including non-breaking spaces) and illegal characters are
    replaced with underscores.
    """
    cleaned_file = f"{out_prefix}.cleaned.fasta"
    with (
        open(input_file, "r", encoding="utf-8") as fin,
        open(cleaned_file, "w", encoding="utf-8") as fout,
    ):
        for line in fin:
            line = line.strip()
            modified_line = "".join(
                "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c for c in line
            )
            fout.write(f"{modified_line}{os.linesep}")
    return cleaned_file


def run_as_subprocess(command, out_file, raise_errors=False):
    """
    Run a command using subprocess.run, piping stdout and stderr to `out_file`.
    """
    with open(out_file, "a", encoding="utf-8") as f:
        result = subprocess.run(
            command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
        )
        if result.returncode != 0:
            command_str = ' '.join(command)
            logging.info(f"\t ERROR: command '{command_str}' failed")
            raise RuntimeError(
                f"subprocess.run of command '{command_str}' encountered error."
                f" Check {out_file} for logs."
            )


def translate_input(input_file, output_prefix, log_file):
    """
    Use `transeq` to get the six-frame translation of the input FASTA. If the input was protein,
    reverse translate it first, overwriting the original input.
    """
    logging.debug("\t...running transeq")

    # Check whether the input file is protein, and should be reverse-translated first
    # We use highly-expressed genes in E. coli, since there's no standard table
    if is_likely_protein(input_file):
        command = ["backtranambig", input_file, input_file]
        run_as_subprocess(command, log_file)

    faa_to_screen = f"{output_prefix}.transeq.faa"
    command = ["transeq", input_file, faa_to_screen, "-frame", "6", "-clean"]
    run_as_subprocess(command, log_file)
    return faa_to_screen


def screen_biorisks(
    input_file, output_prefix, screen_file, log_file, scripts_dir, screen_dbs
):
    """
    Call hmmscan` and `check_biorisk.py` to add biorisk results to `screen_file`.
    """
    logging.debug("\t...running hmmscan")
    biorisk_out = f"{output_prefix}.biorisk.hmmscan"
    command = ["hmmscan", "--domtblout", biorisk_out, screen_dbs.biorisk_db, input_file]
    run_as_subprocess(command, log_file)

    logging.debug("\t...checking hmmscan results")
    command = [
        "python3",
        f"{scripts_dir}/check_biorisk.py",
        "-i",
        biorisk_out,
        "--database",
        screen_dbs.biorisk_dir,
    ]
    run_as_subprocess(command, screen_file)


def screen_proteins(
    input_file,
    out_prefix,
    screen_file,
    tmp_log,
    scripts_dir,
    screen_dbs,
    search_tool,
    threads,
):
    """
    Call `run_blastx.sh` or `run_diamond.sh` followed by `check_reg_path.py` to add regulated
    pathogen protein screening results to `screen_file`.
    """
    logging.debug("\t...running %s", search_tool)
    protein_out = f"{out_prefix}.nr"

    if search_tool == "blastx":
        search_output = f"{out_prefix}.nr.blastx"
        command = [
            f"{scripts_dir}/run_blastx.sh",
            "-q",
            input_file,
            "-d",
            screen_dbs.nr_db,
            "-o",
            protein_out,
            "-t",
            str(threads),
        ]
        run_as_subprocess(command, tmp_log)
    else:  # search_tool == 'diamond':
        search_output = f"{out_prefix}.nr.dmnd"
        processes = 6
        command = [
            f"{scripts_dir}/run_diamond.sh",
            "-i",
            input_file,
            "-d",
            screen_dbs.nr_dir,
            "-o",
            protein_out,
            "-t",
            str(int(threads / processes)),
            "-p",
            str(processes),
        ]
        run_as_subprocess(command, tmp_log)

    if not os.path.exists(search_output):
        raise RuntimeError(f"Protein search failed and {search_output} was not created. Aborting.")

    logging.debug("\t...checking %s results", search_tool)
    # Delete any previous check_reg_path results
    reg_path_coords = f"{out_prefix}.reg_path_coords.csv"
    if os.path.isfile(reg_path_coords):
        os.remove(reg_path_coords)

    command = [
        "python",
        f"{scripts_dir}/check_reg_path.py",
        "-i",
        search_output,
        "-d",
        screen_dbs.dir,
        "-t",
        str(threads),
    ]
    run_as_subprocess(command, screen_file)


def screen_nucleotides(
    input_file, out_prefix, screen_file, scripts_dir, screen_dbs, search_tool, threads
):
    """
    Call `fetch_nc_bits.py`, search noncoding regions with `blastn` and then `check_reg_path.py` to
    screen regulated pathogen nucleotides in noncoding regions (i.e. that would not be found with
    protein search).
    """
    nr_out = out_prefix + ".nr.blastx" if search_tool == "blastx" else "nr.dmnd"
    command = ["python", f"{scripts_dir}/fetch_nc_bits.py", nr_out, input_file]
    run_as_subprocess(command, screen_file)

    # Only screen nucleotides in noncoding regions
    noncoding_fasta = f"{out_prefix}.noncoding.fasta"
    if not os.path.isfile(noncoding_fasta):
        logging.debug(
            "\t...skipping nucleotide search since no noncoding regions fetched"
        )
        return

    # Only run new blastn search if there are no previous results
    nt_output = f"{out_prefix}.nt.blastn"
    if not os.path.isfile(nt_output):
        command = [
            "blastn",
            "-query",
            noncoding_fasta,
            "-db",
            screen_dbs.nt_db,
            "-out",
            nt_output,
            "-outfmt",
            "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send",
            "-max_target_seqs",
            "50",
            "-num_threads",
            "8",
            "-culling_limit",
            "5",
            "-evalue",
            "10",
        ]
        run_as_subprocess(command, screen_file)

    if not os.path.exists(nt_output):
        raise RuntimeError(f"Nucleotide search failed and {nt_output} was not created. Aborting.")

    logging.debug("\t...checking blastn results")
    command = [
        "python",
        f"{scripts_dir}/check_reg_path.py",
        "-i",
        nt_output,
        "-d",
        screen_dbs.dir,
        "-t",
        str(threads),
    ]
    run_as_subprocess(command, screen_file)


def screen_benign(
    input_fasta, input_faa, out_prefix, screen_file, tmp_log, scripts_dir, benign_dir
):
    """
    Call `hmmscan`, `blastn`, and `cmscan` and then pass results to `check_benign.py` to identify
    regions that can be cleared.
    """
    logging.debug("\t...running benign hmmscan")
    hmmscan_out = f"{out_prefix}.benign.hmmscan"
    command = [
        "hmmscan",
        "--domtblout",
        hmmscan_out,
        f"{benign_dir}/benign.hmm",
        input_faa,
    ]
    run_as_subprocess(command, tmp_log)

    logging.debug("\t...running benign blastn")
    command = [
        "blastn",
        "-query",
        input_fasta,
        "-db",
        f"{benign_dir}/benign.fasta",
        "-out",
        f"{out_prefix}.benign.blastn",
        "-outfmt",
        "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send",
        "-evalue",
        "1e-5"
    ]
    run_as_subprocess(command, tmp_log)

    logging.debug("\t...running benign cmscan")
    cmscan_out = f"{out_prefix}.benign.cmscan"
    command = ["cmscan", "--tblout", cmscan_out, f"{benign_dir}/benign.cm", input_fasta]
    run_as_subprocess(command, tmp_log)

    logging.debug("\t...checking benign scan results")
    command = [
        "python3",
        f"{scripts_dir}/check_benign.py",
        "-i",
        out_prefix,
        "--database",
        benign_dir,
    ]
    run_as_subprocess(command, screen_file)

def run(args):
    """
    Wrapper so that args be parsed in main() or commec.py interface.
    """
    input_file = args.fasta_file
    output_prefix = get_output_prefix(input_file, args.output_prefix)
    screen_file = f"{output_prefix}.screen"

    if os.path.exists(screen_file):
        raise RuntimeError(f"Screen output {screen_file} already exists. Aborting.")

    tmp_log = f"{output_prefix}.log.tmp"

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(screen_file, "a")],
    )
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[logging.FileHandler(tmp_log, "a")],
    )

    # Set run parameters from args
    run_parameters = ScreenParameters(
        threads=args.threads,
        protein_search_tool=args.protein_search_tool,
        in_fast_mode=args.fast_mode,
        skip_nt_search=args.skip_nt_search,
        do_cleanup=args.cleanup
    )

    # Check that databases exist
    scripts_dir = os.path.dirname(__file__)
    dbs = ScreenDatabases(args.database_dir, run_parameters.protein_search_tool)
    dbs.validate(run_parameters.in_fast_mode, run_parameters.skip_nt_search)

    # Add the input contents to the log
    shutil.copyfile(input_file, tmp_log)
    fasta_to_screen = get_cleaned_fasta(input_file, output_prefix)

    # Biorisk screen
    logging.info(">> STEP 1: Checking for biorisk genes...")
    faa_to_screen = translate_input(fasta_to_screen, output_prefix, tmp_log)
    screen_biorisks(faa_to_screen, output_prefix, screen_file, tmp_log, scripts_dir, dbs)
    logging.info(
        " STEP 1 completed at %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )

    # Taxonomy screen
    if run_parameters.in_fast_mode:
        logging.info(" >> FAST MODE: Skipping steps 2-3")
    else:
        logging.info(" >> STEP 2: Checking regulated pathogen proteins...")
        screen_proteins(
            fasta_to_screen,
            output_prefix,
            screen_file,
            tmp_log,
            scripts_dir,
            dbs,
            run_parameters.protein_search_tool,
            run_parameters.threads,
        )
        logging.info(
            " STEP 2 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

    if run_parameters.in_fast_mode or run_parameters.skip_nt_search:
        logging.info(" SKIPPING STEP 3: Nucleotide search")
    else:
        logging.info(" >> STEP 3: Checking regulated pathogen nucleotides...")
        screen_nucleotides(
            fasta_to_screen,
            output_prefix,
            screen_file,
            scripts_dir,
            dbs,
            run_parameters.protein_search_tool,
            run_parameters.threads,
        )
        logging.info(
            " STEP 3 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

    # Benign screen
    logging.info(">> STEP 4: Checking any pathogen regions for benign components...")
    screen_benign(
        fasta_to_screen,
        faa_to_screen,
        output_prefix,
        screen_file,
        tmp_log,
        scripts_dir,
        dbs.benign_dir,
    )
    logging.info(
        ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )

    if run_parameters.do_cleanup:
        for pattern in [
            "reg_path_coords.csv",
            "*hmmscan",
            "*blastn",
            "faa",
            "*blastx",
            "*dmnd",
            "*.tmp",
        ]:
            for file in glob.glob(f"{output_prefix}.{pattern}"):
                if os.path.isfile(file):
                    os.remove(file)

def main():
    """
    Main function. Passes args to `run`.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    run(parser.parse_args())

if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        print(f"Runtime error: {e}", file=sys.stderr)
        sys.exit(1)
