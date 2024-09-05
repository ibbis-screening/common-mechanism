#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
Contains information pertinant to commec screen, including input configuration,
and input fasta / query validation and parsing.
"""
import os
import glob
import argparse
from dataclasses import dataclass

from commec.config.query import Query

@dataclass
class ScreenConfiguration:
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
    force_overwrite : bool = False

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
        self.config : ScreenConfiguration = ScreenConfiguration(
            args.threads,
            args.protein_search_tool,
            args.fast_mode,
            args.skip_nt_search,
            args.cleanup,
            args.jobs,
        )

        # TODO: Think about whether logs belong in here, or externally.

        # Outputs
        self.output_prefix = self.get_output_prefix(args.fasta_file, args.output_prefix)
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.tmp_log = f"{self.output_prefix}.log.tmp"

        #Query
        self.query : Query = Query(args.fasta_file)

        # Should this instead be elsewhere... 
        self.db_dir = args.database_dir


    def validate(self) -> bool:
        """
        Make sure all the needed databases exist.
        """
        if os.path.exists(self.output_screen_file):
            raise RuntimeError(f"Screen output {self.output_screen_file} already exists. Aborting.")

        self.query.validate(self.output_prefix)
        return True

    @staticmethod
    def get_output_prefix(input_file, prefix_arg=""):
        """
        Returns a prefix that can be used for all output files.
        
        - If no prefix was given, use the input filename.
        - If a directory was given, use the input filename as file prefix within that directory.
        """
        if not prefix_arg:
            return os.path.splitext(input_file)[0]
        if os.path.isdir(prefix_arg) or prefix_arg.endswith(os.path.sep) or prefix_arg in {".", "..", "~"}:
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
        #return
        if self.config.do_cleanup:
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
        return not self.config.in_fast_mode

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.config.in_fast_mode or self.config.skip_nt_search)

    @property
    def should_do_benign_screening(self) -> bool:
        return True