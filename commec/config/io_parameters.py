#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `ScreenIOParameters` class and associated dataclasses.
Objects responsible for parsing and interpreting user input for
the screen workflow of commec.
"""
import os
import sys
import glob
import argparse
import logging
from dataclasses import dataclass
from typing import Optional
import multiprocessing

import yaml
from yaml.parser import ParserError

from commec.config.query import Query


@dataclass
class ScreenConfig:
    """
    Namespace for optional input parameters for screening; provided by parsing command line
    arguments. Default values, where applicable, are stored here.
    """

    threads: int = 1
    protein_search_tool: str = "blastx"
    in_fast_mode: bool = False
    skip_nt_search: bool = False
    do_cleanup: bool = False
    diamond_jobs: Optional[int] = None
    config_yaml_file: str | os.PathLike = "commec-config.yaml"
    force: bool = False
    resume: bool = False


class ScreenIOParameters:
    """
    Container for input settings constructed from arguments to `screen`.
    """

    def __init__(self, args: argparse.ArgumentParser):
        # Inputs
        self.config: ScreenConfig = ScreenConfig(
            threads=args.threads,
            protein_search_tool=args.protein_search_tool,
            in_fast_mode=args.fast_mode,
            skip_nt_search=args.skip_nt_search,
            do_cleanup=args.cleanup,
            diamond_jobs=args.diamond_jobs,
            config_yaml_file=args.config_yaml.strip(),
            force=args.force,
            resume=args.resume,
        )

        # Outputs
        self.output_prefix = self.get_output_prefix(args.fasta_file, args.output_prefix)
        self.output_screen_file = f"{self.output_prefix}.screen"
        self.tmp_log = f"{self.output_prefix}.log.tmp"

        # Query
        self.query: Query = Query(args.fasta_file)

        # Parse user input paths to local databases:
        self.db_dir = args.database_dir
        self.yaml_configuration = {}
        if os.path.exists(self.config.config_yaml_file):
            self.get_configurations_from_yaml(self.config.config_yaml_file, self.db_dir)
        else:
            print("File not found: " + self.config.config_yaml_file)
            raise FileNotFoundError(
                "No configuration yaml was found at "
                + self.config.config_yaml_file
                + " "
                "if you are using a custom config file, check the path is correct"
            )

        # Check whether a .screen output file already exists.
        if os.path.exists(self.output_screen_file) and not (
            self.config.force or self.config.resume
        ):
            print(
                f"Screen output {self.output_screen_file} already exists. \n"
                "Either use a different output location, or use --force or --resume to override. "
                "\nAborting Screen."
            )
            sys.exit(1)

    def setup(self) -> bool:
        """
        Post-instantiation additonal setup. i.e. setup require output logs.
        """
        # Sanity checks on thread input.
        if self.config.threads > multiprocessing.cpu_count():
            logging.info(
                "Requested allocated threads [%i] is greater"
                " than the detected CPU count of the hardware[%i].",
                self.config.threads,
                multiprocessing.cpu_count(),
            )
        if self.config.threads < 1:
            raise RuntimeError("Number of allocated threads must be at least 1!")

        if (
            self.config.diamond_jobs is not None
            and self.config.protein_search_tool == "blastx"
        ):
            logging.info(
                "WARNING: --jobs is a diamond only parameter! "
                "Specifying -j (--jobs) without also specifying "
                "-p (--protein-search-tool), the protein search "
                'tool as "diamond" will have no effect!'
            )

        self.query.setup(self.output_prefix)
        return True

    def get_configurations_from_yaml(
        self, config_filepath: str, db_dir_override: str = None
    ):
        """
        Read the contents of a YAML file, to see
        if it contains information to override configuration.
        TODO: Candidate for future yaml/json io refactor in future.
        """
        config_from_yaml = None
        try:
            with open(config_filepath, "r", encoding="utf-8") as file:
                config_from_yaml = yaml.safe_load(file)
        except ParserError as e:
            print(
                f"A configuration.yaml file was found ({config_filepath}) "
                "but was invalid as a yaml file:\n",
                e,
            )
            return {}

        # Extract base paths for substitution
        base_paths = {}
        try:
            base_paths = config_from_yaml["base_paths"]
            if db_dir_override:
                base_paths["default"] = db_dir_override

            # Function to recursively replace placeholders
            def recursive_format(d, base_paths):
                if isinstance(d, dict):
                    return {k: recursive_format(v, base_paths) for k, v in d.items()}
                if isinstance(d, str):
                    return d.format(**base_paths)
                return d

            config_from_yaml = recursive_format(config_from_yaml, base_paths)
        except TypeError:
            pass

        self.yaml_configuration = config_from_yaml

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
    def should_do_protein_screening(self) -> bool:
        return not self.config.in_fast_mode

    @property
    def should_do_nucleotide_screening(self) -> bool:
        return not (self.config.in_fast_mode or self.config.skip_nt_search)

    @property
    def should_do_benign_screening(self) -> bool:
        return True
