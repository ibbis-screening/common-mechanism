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
import sys
import pandas as pd

from commec.file_tools import FileTools
from commec.io_parameters import ScreenIOParameters, ScreenInputParameters
from commec.screen_databases import CommecDatabases

from commec.check_biorisk import check_biorisk
from commec.check_benign import check_for_benign
from commec.check_reg_path import check_for_regulated_pathogens
from commec.fetch_nc_bits import fetch_noncoding_regions

from commec.json_io import (
    ScreenData,
    encode_screen_data_to_json,
    #get_screen_data_from_json
)

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
        "-t", 
        "--threads", 
        dest="threads", 
        type=int, 
        default=default_params.threads, 
        help="Threads available"
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

class Screen:
    """ 
    Handles the parsing of input arguments, the control of databases, and 
    the logical flow of the screening process for commec.
    """

    def __init__(self):
        self.params : ScreenIOParameters = None
        self.databases : CommecDatabases = None
        self.scripts_dir : str = os.path.dirname(__file__) # Legacy.
        self.output_screen_data : ScreenData = ScreenData()

    def setup(self, args : argparse.ArgumentParser):
        """ Instantiates and validates parameters, and databases, ready for a run."""
        self.params : ScreenIOParameters = ScreenIOParameters(args)

        #Use print so as not to overwrite an existing .screen file.
        print(" Validating Inputs...")
        self.params.validate()

        # Update filepath.
        #self.output_screen_data_filepath : str =f"{self.params.output_prefix}.results.json"

        logging.basicConfig(level=logging.DEBUG)


        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
            handlers=[logging.StreamHandler(), logging.FileHandler(self.params.output_screen_file, "a")],
        )
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(message)s",
            handlers=[logging.FileHandler(self.params.tmp_log, "a")],
        )

        logging.info(" Validating Inputs...")
        self.databases : CommecDatabases = CommecDatabases(self.params)

        # Add the input contents to the log
        shutil.copyfile(self.params.query.fasta_filepath, self.params.tmp_log)

        # Update screen data output with the Query information.
        self.output_screen_data.query.name = self.params.query.query_description
        self.output_screen_data.query.length = len(self.params.query.aa_raw)
        self.output_screen_data.query.sequence = self.params.query.aa_raw

        # Update screen data output with the commec run information.
        if self.params.should_do_biorisk_screening and not self.databases.biorisk_db is None:
            biorisk_v_info = self.databases.biorisk_db.get_version_information()
            self.output_screen_data.commec_info.biorisk_database_info = biorisk_v_info.version_date

        if self.params.should_do_protein_screening and not self.databases.protein_db is None:
            protein_v_info = self.databases.protein_db.get_version_information()
            self.output_screen_data.commec_info.protein_database_info = protein_v_info.version_date

        if self.params.should_do_nucleotide_screening and not self.databases.nucleotide_db is None:
            nucleotide_v_info = self.databases.nucleotide_db.get_version_information()
            self.output_screen_data.commec_info.nucleotide_database_info = nucleotide_v_info.version_date

        if self.params.should_do_benign_screening and not self.databases.benign_hmm is None:
            benign_v_info = self.databases.benign_hmm.get_version_information()
            self.output_screen_data.commec_info.benign_database_info = benign_v_info.version_date

        encode_screen_data_to_json(self.output_screen_data, self.params.output_json_file)

    def run(self, args : argparse.ArgumentParser):
        """
        Wrapper so that args be parsed in main() or commec.py interface.
        """
        # Perform setup steps.
        self.setup(args)

        # Biorisk screen
        if self.params.should_do_biorisk_screening:
            logging.info(">> STEP 1: Checking for biorisk genes...")
            self.screen_biorisks()
            logging.info(" STEP 1 completed at %s",
                        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        else:
            logging.info(" SKIPPING STEP 1: Biorisk search")

        # Taxonomy screen
        if self.params.should_do_protein_screening:
            logging.info(" >> STEP 2: Checking regulated pathogen proteins...")
            self.screen_proteins()
            #    params.query.fasta_aa_filepath,
            #    params.output_prefix,
            #    params.output_screen_file,
            #    params.tmp_log,
            #    scripts_dir,
            #    params.nr_dir,
            #    params.inputs.search_tool,
            #    params.inputs.threads,
            #    params.inputs.diamond_jobs
            #)
            logging.info(
                " STEP 2 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 2: Protein search")

        if self.params.should_do_nucleotide_screening:
            logging.info(" >> STEP 3: Checking regulated pathogen nucleotides...")
            self.screen_nucleotides()
            #    params.query.cleaned_fasta_filepath,
            #    params.output_prefix,
            #    params.output_screen_file,
            #    scripts_dir,
            #    params,
            #    params.inputs.search_tool,
            #    params.inputs.threads,
            #)
            logging.info(
                " STEP 3 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 3: Nucleotide search")

        # Benign Screen
        if self.params.should_do_benign_screening:
            logging.info(">> STEP 4: Checking any pathogen regions for benign components...")
            self.screen_benign()
            #    params.query.cleaned_fasta_filepath,
            #    params.query.fasta_aa_filepath,
            #    params.output_prefix,
            #    params.output_screen_file,
            #    params.tmp_log,
            #    scripts_dir,
            #    dbs.benign_hmm.db_directory,
            #)
            logging.info(
                ">> STEP 4 completed at %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )
        else:
            logging.info(" SKIPPING STEP 4: Benign search")

        logging.info(
            ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
        self.params.clean()

    # Would like to get rid of this someday, as its currently only used as temporary work around for check_reg_path, which needs work.
    #def run_as_subprocess(self, command, out_file, raise_errors=False):
    #    """
    #    Run a command using subprocess.run, piping stdout and stderr to `out_file`.
    #    """
    #    with open(out_file, "a", encoding="utf-8") as f:
    #        result = subprocess.run(
    #            command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
    #        )
    #        if result.returncode != 0:
    #            command_str = ' '.join(command)
    #            logging.info("\t ERROR: command %s failed", command_str)
    #            raise RuntimeError(
    #                f"subprocess.run of command '{command_str}' encountered error."
    #                f" Check {out_file} for logs."
    #            )

    def screen_biorisks(self):
        """
        Call hmmscan` and `check_biorisk.py` to add biorisk results to `screen_file`.
        """
        logging.debug("\t...running hmmscan")
        self.databases.biorisk_db.screen()
        logging.debug("\t...checking hmmscan results")
        check_biorisk(self.databases.biorisk_db.out_file,
                      self.databases.biorisk_db.db_directory,
                      self.params.output_json_file)

    def screen_proteins(self):
        """
        Call `run_blastx.sh` or `run_diamond.sh` followed by `check_reg_path.py` to add regulated
        pathogen protein screening results to `screen_file`.
        """
        logging.debug("\t...running %s", self.params.inputs.search_tool)
        #protein_out = f"{out_prefix}.nr"

        self.databases.protein_db.screen()

        #if search_tool == "blastx":
        #    search_output = f"{out_prefix}.nr.blastx"
        #    command = [
        #        f"{scripts_dir}/run_blastx.sh",
        #        "-q",
        #        input_file,
        #        "-d",
        #        screen_dbs.nr_db,
        #        "-o",
        #        protein_out,
        #        "-t",
        #        str(threads),
        #    ]
        #    run_as_subprocess(command, tmp_log)
        #else:  # search_tool == 'diamond':
        #    search_output = f"{out_prefix}.nr.dmnd"
        #    command = [
        #        f"{scripts_dir}/run_diamond.sh",
        #        "-i",
        #        input_file,
        #        "-d",
        #        screen_dbs.nr_dir,
        #        "-o",
        #        protein_out,
        #        "-t",
        #        str(threads)
        #    ]
        #    if jobs is not None:
        #        command.extend(['-j', str(jobs)])
        #    run_as_subprocess(command, tmp_log)
#
        if not self.databases.protein_db.check_output(): #os.path.exists(search_output):
            raise RuntimeError(f"Protein search failed and {self.databases.protein_db.out_file} was not created. Aborting.")

        logging.debug("\t...checking %s results", self.params.inputs.search_tool)
        # Delete any previous check_reg_path results
        reg_path_coords = f"{self.params.output_prefix}.reg_path_coords.csv"
        if os.path.isfile(reg_path_coords):
            os.remove(reg_path_coords)

        #command = [
        #    "python",
        #    f"{self.scripts_dir}/check_reg_path.py",
        #    "-i",
        #    self.databases.protein_db.out_file,
        #    "-d",
        #    self.params.db_dir,
        #    "-t",
        #    str(self.params.inputs.threads),
        #]
        #self.run_as_subprocess(command, self.params.output_screen_file)
        check_for_regulated_pathogens(self.databases.protein_db.out_file,
                                       self.params.db_dir,
                                       str(self.params.inputs.threads),
                                       self.params.output_json_file)


    def screen_nucleotides(self
        #input_file, out_prefix, screen_file, scripts_dir, screen_dbs, search_tool, threads
    ):
        """
        Call `fetch_nc_bits.py`, search noncoding regions with `blastn` and then `check_reg_path.py` to
        screen regulated pathogen nucleotides in noncoding regions (i.e. that would not be found with
        protein search).
        """
        #nr_out = out_prefix + ".nr.blastx" if search_tool == "blastx" else "nr.dmnd"
        #command = ["python", f"{scripts_dir}/fetch_nc_bits.py", nr_out, input_file]
        #run_as_subprocess(command, screen_file)

        fetch_noncoding_regions(self.databases.protein_db.out_file, 
                                self.params.query.cleaned_fasta_filepath)

        # Only screen nucleotides in noncoding regions
        noncoding_fasta = f"{self.params.output_prefix}.noncoding.fasta"
        if not os.path.isfile(noncoding_fasta):
            logging.debug(
                "\t...skipping nucleotide search since no noncoding regions fetched"
            )
            return

        # Only run new blastn search if there are no previous results
        if not self.databases.nucleotide_db.check_output():
            self.databases.nucleotide_db.screen()

        #nt_output = f"{out_prefix}.nt.blastn"
        #if not os.path.isfile(nt_output):
        #    command = [
        #        "blastn",
        #        "-query",
        #        noncoding_fasta,
        #        "-db",
        #        screen_dbs.nt_db,
        #        "-out",
        #        nt_output,
        #        "-outfmt",
        #        "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send",
        #        "-max_target_seqs",
        #        "50",
        #        "-num_threads",
        #        "8",
        #        "-culling_limit",
        #        "5",
        #        "-evalue",
        #        "10",
        #    ]
        #    run_as_subprocess(command, screen_file)

        if not self.databases.nucleotide_db.check_output(): #os.path.exists(nt_output):
            raise RuntimeError(f"Nucleotide search failed and {self.databases.nucleotide_db.out_file} was not created. Aborting.")

        logging.debug("\t...checking blastn results")
        #command = [
        #    "python",
        #    f"{self.scripts_dir}/check_reg_path.py",
        #    "-i",
        #    self.databases.nucleotide_db.out_file,
        #    "-d",
        #    self.params.db_dir,
        #    "-t",
        #    str(self.params.inputs.threads),
        #]
        #self.run_as_subprocess(command, self.params.output_screen_file)

        check_for_regulated_pathogens(self.databases.nucleotide_db.out_file,
                                       self.params.db_dir,
                                       str(self.params.inputs.threads),
                                       self.params.output_json_file)


    def screen_benign(self#,
        #input_fasta, input_faa, out_prefix, screen_file, tmp_log, scripts_dir, benign_dir
    ):
        """
        Call `hmmscan`, `blastn`, and `cmscan` and then pass results to `check_benign.py` to identify
        regions that can be cleared.
        """
        logging.debug("\t...running benign hmmscan")
        self.databases.benign_hmm.screen()
        logging.debug("\t...running benign blastn")
        self.databases.benign_blastn.screen()
        logging.debug("\t...running benign cmscan")
        self.databases.benign_cmscan.screen()

        sample_name = self.params.output_prefix # Note currently check_for_benign hard codes .benign.hmmscan onto this.
        if not os.path.exists(sample_name + ".reg_path_coords.csv"):
            logging.info("\t...no regulated regions to clear\n")
            return
        
        coords = pd.read_csv(sample_name + ".reg_path_coords.csv")
        benign_desc =  pd.read_csv(self.databases.benign_hmm.db_directory + "/benign_annotations.tsv", sep="\t")
        
        logging.debug("\t...checking benign scan results")
        check_for_benign(sample_name, coords, benign_desc, self.params.output_json_file)

        #logging.debug("\t...running benign hmmscan")
        #hmmscan_out = f"{out_prefix}.benign.hmmscan"
        #command = [
        #    "hmmscan",
        #    "--domtblout",
        #    hmmscan_out,
        #    f"{benign_dir}/benign.hmm",
        #    input_faa,
        #]
        #run_as_subprocess(command, tmp_log)
#
        #logging.debug("\t...running benign blastn")
        #command = [
        #    "blastn",
        #    "-query",
        #    input_fasta,
        #    "-db",
        #    f"{benign_dir}/benign.fasta",
        #    "-out",
        #    f"{out_prefix}.benign.blastn",
        #    "-outfmt",
        #    "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send",
        #    "-evalue",
        #    "1e-5"
        #]
        #run_as_subprocess(command, tmp_log)
#
        #logging.debug("\t...running benign cmscan")
        #cmscan_out = f"{out_prefix}.benign.cmscan"
        #command = ["cmscan", "--tblout", cmscan_out, f"{benign_dir}/benign.cm", input_fasta]
        #run_as_subprocess(command, tmp_log)

        #command = [
        #    "python3",
        #    f"{scripts_dir}/check_benign.py",
        #    "-i",
        #    self.params.output_prefix,
        #    "--database",
        #    benign_dir,
        #]
        #run_as_subprocess(command, screen_file)

def run(args : argparse.ArgumentParser):
    """
    Entry point from commec. Passes args to Screen object, and runs.
    """
    my_screen : Screen = Screen()
    my_screen.run(args)

def main():
    """
    Main function. Passes args to Screen object, which then runs.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        print(f"Runtime error: {e}", file=sys.stderr)
        sys.exit(1)
