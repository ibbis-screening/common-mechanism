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
import logging
import os
import shutil
import subprocess
import sys

from commec.file_tools import FileTools
from commec.io_parameters import ScreenIOParameters, ScreenInputParameters
from commec.screen_databases import CommecDatabases, DatabaseHandler

DESCRIPTION = "Run Common Mechanism screening on an input FASTA."

def add_args(parser : argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """

    default_params : ScreenInputParameters = ScreenInputParameters()

    parser.add_argument(
        dest="fasta_file",
        type=FileTools.file_arg,
        help="FASTA file to screen"
    )
    parser.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=FileTools.directory_arg,
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
        "-t", "--threads", dest="threads", type=int, default=default_params.threads, help="Threads available"
    )
    parser.add_argument(
        "-j",
        "--jobs",
        dest="jobs",
        type=int,
        help="number of diamond runs to do in parallel (optional, defaults to # CPUs / THREADS)"
    )
    parser.add_argument(
        "-p",
        "--protein-search-tool",
        dest="protein_search_tool",
        choices=["blastx", "diamond"],
        default=default_params.search_tool,
        help="Tool for homology search to identify regulated pathogen proteins",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        dest="jobs",
        type=int,
        help="number of diamond runs to do in parallel (optional, defaults to # CPUs / THREADS)"
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
            logging.info("\t ERROR: command %s failed", command_str)
            raise RuntimeError(
                f"subprocess.run of command '{command_str}' encountered error."
                f" Check {out_file} for logs."
            )

def screen_biorisks(params: ScreenIOParameters, database : DatabaseHandler):
    """
    Call hmmscan` and `check_biorisk.py` to add biorisk results to `screen_file`.
    """
    logging.debug("\t...running hmmscan")

    #biorisk_in = params.query.fasta_aa_filepath
    #biorisk_out = f"{params.output_prefix}.biorisk.hmmscan"
    #biorisk_dir = os.path.join(params.db_dir, "biorisk_db")
    #biorisk_db = os.path.join(biorisk_dir, "biorisk.hmm")

    #my_biorisk_screener = HMMDataBase(biorisk_dir, biorisk_db, biorisk_in, biorisk_out)
    database.screen()

    logging.debug("\t...checking hmmscan results")
    command = [
        "python3",
        f"{os.path.dirname(__file__)}/check_biorisk.py",
        "-i",
        database.out_file,
        "--database",
        database.db_directory,
    ]
    run_as_subprocess(command, params.output_screen_file)


def screen_proteins(
    input_file,
    out_prefix,
    screen_file,
    tmp_log,
    scripts_dir,
    screen_dbs,
    search_tool,
    threads,
    jobs=None
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
        command = [
            f"{scripts_dir}/run_diamond.sh",
            "-i",
            input_file,
            "-d",
            screen_dbs.nr_dir,
            "-o",
            protein_out,
            "-t",
            str(threads)
        ]
        if jobs is not None:
            command.extend(['-j', str(jobs)])
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
    params : ScreenIOParameters = ScreenIOParameters(args)
    #Use print so as not to overwrite an existing .screen file.
    print(" Validating Inputs...") 
    params.validate()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(params.output_screen_file, "a")],
    )

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[logging.FileHandler(params.tmp_log, "a")],
    )

    # Check that databases exist
    scripts_dir = os.path.dirname(__file__)

    dbs : CommecDatabases = CommecDatabases(params)

    # Add the input contents to the log
    shutil.copyfile(params.query.fasta_filepath, params.tmp_log)

    # Biorisk screen
    if params.should_do_biorisk_screening:
        logging.info(">> STEP 1: Checking for biorisk genes...")
        screen_biorisks(params, dbs.biorisk_db)
        logging.info(" STEP 1 completed at %s",
                     datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    else:
        logging.info(" SKIPPING STEP 1: Biorisk search")

    # Taxonomy screen
    if params.should_do_protein_screening:
        logging.info(" >> STEP 2: Checking regulated pathogen proteins...")
        screen_proteins(
            params.query.fasta_aa_filepath,
            params.output_prefix,
            params.output_screen_file,
            params.tmp_log,
            scripts_dir,
            params.nr_dir,
            params.inputs.search_tool,
            params.inputs.threads,
            params.inputs.diamond_jobs
        )
        logging.info(
            " STEP 2 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
    else:
        logging.info(" SKIPPING STEP 2: Protein search")

    if params.should_do_nucleotide_screening:
        logging.info(" >> STEP 3: Checking regulated pathogen nucleotides...")
        screen_nucleotides(
            params.query.cleaned_fasta_filepath,
            params.output_prefix,
            params.output_screen_file,
            scripts_dir,
            params,
            params.inputs.search_tool,
            params.inputs.threads,
        )
        logging.info(
            " STEP 3 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
    else:
        logging.info(" SKIPPING STEP 3: Nucleotide search")

    # Benign Screen
    if params.should_do_benign_screening:
        logging.info(">> STEP 4: Checking any pathogen regions for benign components...")
        screen_benign(
            params.query.cleaned_fasta_filepath,
            params.query.fasta_aa_filepath,
            params.output_prefix,
            params.output_screen_file,
            params.tmp_log,
            scripts_dir,
            params.benign_dir,
        )
        logging.info(
            ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
    else:
        logging.info(" SKIPPING STEP 4: Benign search")

    params.clean()

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
