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
from commec.databases.screen_databases import CommecDatabases

from commec.check_biorisk import check_biorisk
from commec.check_benign import check_for_benign
from commec.check_reg_path import check_for_regulated_pathogens
from commec.fetch_nc_bits import fetch_noncoding_regions

from commec.json_io import (
    ScreenData,
    encode_screen_data_to_json,
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

        print(self.params.output_screen_file)

        logging.basicConfig(level=logging.DEBUG,
                            format="%(message)s",
                            handlers=[logging.StreamHandler(), logging.FileHandler(self.params.output_screen_file, "a")],
                            )


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
        # TODO: Parse the input query, separate out the input names, and 
        #self.output_screen_data.query.name = self.params.query.query_description
        #self.output_screen_data.query.length = len(self.params.query.aa_raw)
        #self.output_screen_data.query.sequence = self.params.query.aa_raw

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

        # Taxonomy screen (Protein)
        if self.params.should_do_protein_screening:
            logging.info(" >> STEP 2: Checking regulated pathogen proteins...")
            self.screen_proteins()
            logging.info(
                " STEP 2 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 2: Protein search")

        # Taxonomy screen (Nucleotide)
        if self.params.should_do_nucleotide_screening:
            logging.info(" >> STEP 3: Checking regulated pathogen nucleotides...")
            self.screen_nucleotides()
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
            logging.info(
                ">> STEP 4 completed at %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )
        else:
            logging.info(" SKIPPING STEP 4: Benign search")

        logging.info(
            ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
        self.params.clean()

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
        self.databases.protein_db.screen()
        if not self.databases.protein_db.check_output(): #os.path.exists(search_output):
            raise RuntimeError(f"Protein search failed and {self.databases.protein_db.out_file} was not created. Aborting.")

        logging.debug("\t...checking %s results", self.params.inputs.search_tool)
        # Delete any previous check_reg_path results
        reg_path_coords = f"{self.params.output_prefix}.reg_path_coords.csv"
        if os.path.isfile(reg_path_coords):
            os.remove(reg_path_coords)

        check_for_regulated_pathogens(self.databases.protein_db.out_file,
                                       self.params.db_dir,
                                       str(self.params.inputs.threads),
                                       self.params.output_json_file)


    def screen_nucleotides(self):
        """
        Call `fetch_nc_bits.py`, search noncoding regions with `blastn` and then `check_reg_path.py` to
        screen regulated pathogen nucleotides in noncoding regions (i.e. that would not be found with
        protein search).
        """
        # Only screen nucleotides in noncoding regions
        fetch_noncoding_regions(self.databases.protein_db.out_file,
                                self.params.query.cleaned_fasta_filepath)
        
        noncoding_fasta = f"{self.params.output_prefix}.noncoding.fasta" # TODO: This should be passed into fetch_noncoding_regions.

        if not os.path.isfile(noncoding_fasta):
            logging.debug(
                "\t...skipping nucleotide search since no noncoding regions fetched"
            )
            return

        # Only run new blastn search if there are no previous results
        if not self.databases.nucleotide_db.check_output():
            self.databases.nucleotide_db.screen()

        if not self.databases.nucleotide_db.check_output():
            raise RuntimeError(f"Nucleotide search failed and {self.databases.nucleotide_db.out_file} was not created. Aborting.")

        logging.debug("\t...checking blastn results")
        check_for_regulated_pathogens(self.databases.nucleotide_db.out_file,
                                       self.params.db_dir,
                                       str(self.params.inputs.threads),
                                       self.params.output_json_file)

    def screen_benign(self):
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

def run(args : argparse.ArgumentParser):
    """
    Entry point from commec main. Passes args to Screen object, and runs.
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
