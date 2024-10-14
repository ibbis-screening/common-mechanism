#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Run Common Mechanism screening on an input FASTA.

Screening involves (up to) four steps:

  1. Biorisk scan:      HMM-based scan for matches to a custom database of biorisk sequences.
  2. Protein search:    protein homology search for best matches to regulated pathogens.
  3. Nucleotide search: nucleotide homology search for best matches to regulated pathogens.
  4. Benign scan:       three different scans (against conserved proteins, housekeeping RNAs, and
                        synbio parts) to see if hits identified in homology search can be cleared. 

In "fast" mode, only the biorisk scan is run. By default, all four steps are run, but the nucleotide
search is only run for regions that do not have any protein hits with a high sequence identity. The
benign search is not permitted to clear biorisk scan hits, only protein or nucleotide hits. Whether
or not a homology scan hit is from a regulated pathogen is determined by referencing the taxonomy
ids assoicated with each accession that returns a hit, then looking at their lineages.

positional arguments:
  fasta_file            FASTA file to screen

options:
  -h, --help            show this help message and exit
  -d DATABASE_DIR, --databases DATABASE_DIR
                        Path to directory containing reference databases (e.g. taxonomy, protein,
                        HMM)
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        Prefix for output files. Can be a string (interpreted as output basename) or
                        a directory (in which case the output file names will be determined from the
                        input FASTA)
  -t THREADS, --threads THREADS
                        Number of CPU threads to use. Passed to search tools.
  -p {blastx,diamond}, --protein-search-tool {blastx,diamond}
                        Tool for homology search to identify regulated pathogen proteins
  -j DIAMOND_JOBS, --diamond-jobs DIAMOND_JOBS
                        Diamond-only: number of runs to do in parallel
  -f, --fast            Run in fast mode and skip homology search
  -n, --skip-nt         Skip nucleotide search, even for regions where no protein hits were found
  -c, --cleanup         Delete intermediate files

"""
import argparse
import datetime
import time
import logging
import os
import shutil
import sys
import pandas as pd

from commec.utils.file_utils import file_arg, directory_arg
from commec.config.io_parameters import ScreenIOParameters, ScreenConfig
from commec.config.screen_tools import ScreenTools

from commec.screeners.check_biorisk import check_biorisk, update_biorisk_data_from_database
from commec.screeners.check_benign import check_for_benign
from commec.screeners.check_reg_path import check_for_regulated_pathogens, update_taxonomic_data_from_database
from commec.tools.fetch_nc_bits import fetch_noncoding_regions

from commec.config.json_io import (
    ScreenData,
    CommecScreenStep,
    QueryData,
    encode_screen_data_to_json
)

DESCRIPTION = "Run Common Mechanism screening on an input FASTA."


def add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """
    default_config: ScreenConfig = ScreenConfig()

    parser.add_argument(dest="fasta_file", type=file_arg, help="FASTA file to screen")
    parser.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        required=True,
        help="Path to directory containing reference databases (e.g. taxonomy, protein, HMM)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_prefix",
        help="Prefix for output files. Can be a string (interpreted as output basename) or a" +
         " directory (in which case the output file names will be determined from the input FASTA)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=default_config.threads,
        help="Number of CPU threads to use. Passed to search tools.",
    )
    parser.add_argument(
        "-p",
        "--protein-search-tool",
        dest="protein_search_tool",
        choices=["blastx", "diamond"],
        default=default_config.protein_search_tool,
        help="Tool for homology search to identify regulated pathogen proteins",
    )
    parser.add_argument(
        "-j",
        "--diamond-jobs",
        dest="diamond_jobs",
        type=int,
        default=default_config.diamond_jobs,
        help="Diamond-only: number of runs to do in parallel",
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
        help="Skip nucleotide search, even for regions where no protein hits were found",
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
        self.database_tools : ScreenTools = None
        self.scripts_dir : str = os.path.dirname(__file__) # Legacy.
        self.screen_data : ScreenData = ScreenData()
        self.start_time = time.time()

    def setup(self, args: argparse.ArgumentParser):
        """Instantiates and validates parameters, and databases, ready for a run."""
        self.params: ScreenIOParameters = ScreenIOParameters(args)
        self.params.setup()

        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(self.params.output_screen_file, "a"),
                logging.FileHandler(self.params.tmp_log, "a"),
            ],
        )
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(message)s",
            handlers=[logging.FileHandler(self.params.tmp_log, "a")],
        )

        logging.info(" Validating Inputs...")
        self.params.setup()
        self.database_tools: ScreenTools = ScreenTools(self.params)
        self.params.query.translate_query()

        # Add the input contents to the log
        shutil.copyfile(self.params.query.input_fasta_path, self.params.tmp_log)

        # Initialise the json file:
        for query in self.params.query.raw:
            self.screen_data.queries.append(QueryData(
                query.name,
                len(str(query.seq)),
                str(query.seq))
                )
            
        self.screen_data.commec_info.biorisk_database_info = self.database_tools.biorisk_hmm.get_version_information()

        if self.params.should_do_protein_screening:
            self.screen_data.commec_info.protein_database_info = self.database_tools.regulated_protein.get_version_information()

        if self.params.should_do_nucleotide_screening:
            self.screen_data.commec_info.nucleotide_database_info = self.database_tools.regulated_nt.get_version_information()

        if self.params.should_do_benign_screening:
            self.screen_data.commec_info.benign_protein_database_info = self.database_tools.benign_hmm.get_version_information()
            self.screen_data.commec_info.benign_rna_database_info = self.database_tools.benign_blastn.get_version_information()
            self.screen_data.commec_info.benign_synbio_database_info = self.database_tools.benign_cmscan.get_version_information()

        self.screen_data.commec_info.date_run = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Less of this, and more just passing the object in future - this might take too long in larger queries.
        #encode_screen_data_to_json(self.screen_data, self.params.output_json)
        

    def run(self, args : argparse.ArgumentParser):
        """
        Wrapper so that args be parsed in main() or commec.py interface.
        """
        # Perform setup steps.
        self.setup(args)

        # Biorisk screen
        logging.info(">> STEP 1: Checking for biorisk genes...")
        self.screen_biorisks()
        logging.info(
            " STEP 1 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

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
            logging.info(
                ">> STEP 4: Checking any pathogen regions for benign components..."
            )
            self.screen_benign()
            logging.info(
                ">> STEP 4 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 4: Benign search")

        logging.info(
            ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )

        #self.screen_data = get_screen_data_from_json(self.params.output_json)
        time_taken = (time.time() - self.start_time)
        # Convert the elapsed time to hours, minutes, and seconds
        hours, rem = divmod(time_taken, 3600)
        minutes, seconds = divmod(rem, 60)
        self.screen_data.commec_info.time_taken = f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"
        self.screen_data.update()
        encode_screen_data_to_json(self.screen_data, self.params.output_json)

        self.params.clean()

    def screen_biorisks(self):
        """
        Call hmmscan` and `check_biorisk.py` to add biorisk results to `screen_file`.
        """
        logging.debug("\t...running hmmscan")
        self.database_tools.biorisk_hmm.search()
        logging.debug("\t...checking hmmscan results")
        check_biorisk(
            self.database_tools.biorisk_hmm.out_file,
            self.database_tools.biorisk_hmm.db_directory
        )
        update_biorisk_data_from_database(self.database_tools.biorisk_hmm, self.screen_data)

    def screen_proteins(self):
        """
        Call `run_blastx.sh` or `run_diamond.sh` followed by `check_reg_path.py` to add regulated
        pathogen protein screening results to `screen_file`.
        """
        logging.debug("\t...running %s", self.params.config.protein_search_tool)
        self.database_tools.regulated_protein.search()
        if (
            not self.database_tools.regulated_protein.check_output()
        ):  # os.path.exists(search_output):
            raise RuntimeError(
                "Protein search failed and "
                f"{self.database_tools.regulated_protein.out_file}"
                " was not created. Aborting."
            )

        logging.debug(
            "\t...checking %s results", self.params.config.protein_search_tool
        )
        # Delete any previous check_reg_path results
        reg_path_coords = f"{self.params.output_prefix}.reg_path_coords.csv"
        if os.path.isfile(reg_path_coords):
            os.remove(reg_path_coords)

        check_for_regulated_pathogens(
            self.database_tools.regulated_protein.out_file,
            self.params.db_dir,
            str(self.params.config.threads),
        )
        update_taxonomic_data_from_database(self.database_tools.regulated_protein, 
                                            self.database_tools.benign_blastn, 
                                            self.database_tools.biorisk_hmm,
                                            self.params.db_dir + "/taxonomy/",
                                            self.screen_data,
                                            CommecScreenStep.TAXONOMY_AA,
                                            self.params.config.threads)


    def screen_nucleotides(self):
        """
        Call `fetch_nc_bits.py`, search noncoding regions with `blastn` and
        then `check_reg_path.py` to screen regulated pathogen nucleotides in
        noncoding regions (i.e. that would not be found with protein search).
        """
        # Only screen nucleotides in noncoding regions
        fetch_noncoding_regions(
            self.database_tools.regulated_protein.out_file, self.params.query.nt_path
        )
        noncoding_fasta = f"{self.params.output_prefix}.noncoding.fasta"

        if not os.path.isfile(noncoding_fasta):
            logging.debug(
                "\t...skipping nucleotide search since no noncoding regions fetched"
            )
            return

        # Only run new blastn search if there are no previous results
        if not self.database_tools.regulated_nt.check_output():
            self.database_tools.regulated_nt.search()

        if not self.database_tools.regulated_nt.check_output():
            raise RuntimeError(
                "Nucleotide search failed and "
                f"{self.database_tools.regulated_nt.out_file}"
                " was not created. Aborting."
            )

        logging.debug("\t...checking blastn results")
        check_for_regulated_pathogens(
            self.database_tools.regulated_nt.out_file,
            self.params.db_dir,
            str(self.params.config.threads),
        )
        update_taxonomic_data_from_database(self.database_tools.regulated_nt, 
                                            self.database_tools.benign_blastn, 
                                            self.database_tools.biorisk_hmm,
                                            self.params.db_dir + "/taxonomy/",
                                            self.screen_data,
                                            CommecScreenStep.TAXONOMY_AA,
                                            self.params.config.threads)

    def screen_benign(self):
        """
        Call `hmmscan`, `blastn`, and `cmscan` and then pass results
        to `check_benign.py` to identify regions that can be cleared.
        """
        sample_name = self.params.output_prefix
        if not os.path.exists(sample_name + ".reg_path_coords.csv"):
            logging.info("\t...no regulated regions to clear\n")
            return

        logging.debug("\t...running benign hmmscan")
        self.database_tools.benign_hmm.search()
        logging.debug("\t...running benign blastn")
        self.database_tools.benign_blastn.search()
        logging.debug("\t...running benign cmscan")
        self.database_tools.benign_cmscan.search()

        coords = pd.read_csv(sample_name + ".reg_path_coords.csv")
        benign_desc = pd.read_csv(
            self.database_tools.benign_hmm.db_directory + "/benign_annotations.tsv",
            sep="\t",
        )

        logging.debug("\t...checking benign scan results")

        # Note currently check_for_benign hard codes .benign.hmmscan,
        # in future parse, and grab from search handler instead.
        check_for_benign(sample_name, coords, benign_desc)


def run(args: argparse.ArgumentParser):
    """
    Entry point from commec main. Passes args to Screen object, and runs.
    """
    my_screen: Screen = Screen()
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
