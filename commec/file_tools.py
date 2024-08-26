#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for FileTools, containing Static functions useful for dealing with common file parsing tasks.
"""

import argparse
import os

class FileTools():
    """
    Static function container, useful for dealing with common file parsing tasks.
    """

    @staticmethod
    def is_empty(filepath: str) -> bool:
        """
        is_empty

        usage: check that a file is empty
        input:
        - name of file
        """
        try:
            abspath = os.path.abspath(os.path.expanduser(filepath))
            return os.path.getsize(abspath) == 0
        except OSError:
            # If there is an error (including FileNotFoundError) consider it empty
            return True

    @staticmethod
    def has_hits(filepath: str) -> bool:
        """
        has_hits
        usage: check to see if the file contains any hits (lines that don't start with #)
        input:
        - path to file
        """
        try:
            with open(filepath, "r", encoding="utf-8") as file:
                for line in file:
                    # Strip leading and trailing whitespace and check the first character
                    if not line.strip().startswith("#"):
                        # Found a hit!
                        return True
            return False
        except FileNotFoundError:
            # The file does not exist
            return False

    @staticmethod
    def directory_arg(path):
        """Raise ArgumentTypeError if `path` is not a directory."""
        if not os.path.isdir(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid directory path")
        return path

    @staticmethod
    def file_arg(path):
        """Raise ArgumentTypeError if `path` is not a file."""
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid file")
        if not os.path.getsize(path) > 0:
            raise argparse.ArgumentTypeError(f"{path} is an empty file")
        return path
    
    @staticmethod
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
