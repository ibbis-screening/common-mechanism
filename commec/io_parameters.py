#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
import glob
import argparse
from dataclasses import dataclass
import subprocess

from commec.file_tools import FileTools

@dataclass
class ScreenInputParameters:
    """
    Namespace for optional input parameters for screening; provided by parsing arguments. 
    All have default values, which are set here, and are reflected into argparse.
    """
    threads: int = 1
    search_tool: str = "blastx"
    in_fast_mode: bool = False
    skip_nt_search: bool = False
    do_cleanup: bool = False
    diamond_jobs : int = None

@dataclass
class Query:
    """
    Container for all information related to a query. Self calculates AA version.
    Potential to expand in future to reverse translate from AA too.
    """
    fasta_filepath : str

    fasta_aa_filepath : str = ""
    cleaned_fasta_filepath : str = ""

    # Unused for now.
    query_description : str = ""
    nt_raw : str = "" # The full query in string format eg: ATGGCTGACGATCGATCGACTCG
    aa_raw : str = "" # The full query in string format eg: MAPTGWKPHGDCGHHHHHHGDPM

    def validate(self, output_prefix : str):
        """ Translate or reverse translate query, so we have it ready in AA or NT format. """
        self.cleaned_fasta_filepath = FileTools.get_cleaned_fasta(self.fasta_filepath, output_prefix)
        self.fasta_aa_filepath = f"{output_prefix}.transeq.faa"
        command = ["transeq", self.cleaned_fasta_filepath, self.fasta_aa_filepath, "-frame", "6", "-clean"]

        # hard coded into subprocess.run for now, as have yet to decide where run_as_subprocess should be.
        # TODO: Update this.
        try:
            result = subprocess.run(command, check=True)
            if result.returncode != 0:
                command_str = ' '.join(command)
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                )
        except RuntimeError as error:
            raise RuntimeError("Input FASTA {fasta_to_screen} could not be translated.") from error

class ScreenIOParameters():
    """
    Container for all input settings, and for all locations and directory structures of inputs / outputs.
    Is constructed with the argsparser outputs, and validates all input settings, and creates secondary
    input settings - such as cleaning fastas, and running transeq.
    """

    def __init__(self, args : argparse.ArgumentParser):
        # Ensure required arguments are present
        required_args = ['fasta_file', 'threads', 'protein_search_tool', 'fast_mode',
                         'skip_nt_search', 'cleanup', 'output_prefix', 'database_dir', 'jobs']
        for arg in required_args:
            if not hasattr(args, arg):
                raise ValueError(f"Missing required argument: {arg}")

        # Inputs
        self.inputs : ScreenInputParameters = ScreenInputParameters(
            args.threads,
            args.protein_search_tool,
            args.fast_mode,
            args.skip_nt_search,
            args.cleanup,
            args.jobs
        )

        # TODO: Think about whether logs belong in here, or externally.

        # Outputs
        self.output_prefix = self.get_output_prefix(args.fasta_file, args.output_prefix)
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.output_json_file = f"{self.output_prefix}.results.json"
        self.tmp_log = f"{self.output_prefix}.log.tmp"

        #Query
        self.query : Query = Query(args.fasta_file, f"{self.output_prefix}.transeq.faa")

        # Should this instead be elsewhere.
        self.db_dir = args.database_dir


    def validate(self) -> bool:
        """
        Make sure all the needed databases exist.
        """
        if os.path.exists(self.output_screen_file):
            raise RuntimeError(f"Screen output {self.output_screen_file} already exists. Aborting.")

        self.query.validate(self.output_prefix)
        return True

    @property
    def nr_dir(self):
        """
        Get the directory containing the needed NR database (which depends on the tool used).
        """
        if self.inputs.search_tool == "blastx":
            return os.path.join(self.db_dir, "nr_blast")
        elif self.inputs.search_tool == "diamond":
            return os.path.join(self.db_dir, "nr_dmnd")
        else:
            raise ValueError(
                "Protein search tool must be either 'blastx' or 'diamond'."
            )

    @property
    def nr_db(self):
        """
        Get the directory of the needed NR database (which depends on the tool used).
        """
        if self.inputs.search_tool == "blastx":
            return os.path.join(self.nr_dir, "nr")
        elif self.inputs.search_tool == "diamond":
            return os.path.join(self.nr_dir, "nr.dmnd")
        else:
            raise ValueError(
                "Protein search tool must be either 'blastx' or 'diamond'."
            )

    @staticmethod
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

    def clean(self):
        """
        Tidy up some directories and erroneous files after a run.
        """
        return
        if self.inputs.do_cleanup:
            for pattern in [
                "reg_path_coords.csv",
                "*hmmscan",
                "*blastn",
                "faa",
                "*blastx",
                "*dmnd",
                "*.tmp",
            ]:
                for file in glob.glob(f"{self.output_prefix}.{pattern}"):
                    if os.path.isfile(file):
                        os.remove(file)

    @property
    def should_do_biorisk_screening(self) -> bool:
        return True

    @property
    def should_do_protein_screening(self) -> bool:
        return not self.inputs.in_fast_mode

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.inputs.in_fast_mode or self.inputs.skip_nt_search)

    @property
    def should_do_benign_screening(self) -> bool:
        return True