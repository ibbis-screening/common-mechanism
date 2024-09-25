#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `ScreenIOParameters` class and associated dataclasses.
"""
import os
import glob
import argparse
from dataclasses import dataclass
from typing import Optional

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
    diamond_jobs: Optional[int] = None


class ScreenIOParameters:
    """
    Container for input settings constructed from arguments to `screen`.
    """
    def __init__(self, args: argparse.ArgumentParser):
        # Inputs
        self.config: ScreenConfiguration = ScreenConfiguration(
            args.threads,
            args.protein_search_tool,
            args.fast_mode,
            args.skip_nt_search,
            args.cleanup,
            args.jobs,
        )

        # Outputs
        self.output_prefix = self.get_output_prefix(args.fasta_file, args.output_prefix)
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.tmp_log = f"{self.output_prefix}.log.tmp"

        # Query
        self.query: Query = Query(args.fasta_file)

        # Storage of user input, used by screen_tools.
        self.db_dir = args.database_dir

    def setup(self) -> bool:
        """
        Post-instantiation additonal setup. Or setup which requires logging.
        """
        #TODO: add -force flag to enable overwriting?
        if os.path.exists(self.output_screen_file):
            raise RuntimeError(
                f"Screen output {self.output_screen_file} already exists. Aborting."
            )

        self.query.setup(self.output_prefix)
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
        if (
            os.path.isdir(prefix_arg)
            or prefix_arg.endswith(os.path.sep)
            or prefix_arg in {".", "..", "~"}
        ):
            os.makedirs(prefix_arg, exist_ok=True)
            # Use the input filename as a prefix within that directory (stripping out the path)
            input_name = os.path.splitext(os.path.basename(input_file))[0]
            return prefix_arg + input_name
        # Existing, non-path prefixes can be used as-is
        return prefix_arg

    def clean(self):
        """
        Tidy up directories and temporary files after a run.
        """
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
