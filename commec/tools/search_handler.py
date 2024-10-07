#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Abstract base class defining a shared interface for search tools.
"""
from abc import ABC, abstractmethod
import os
from dataclasses import dataclass
import subprocess
import logging


@dataclass
class SearchToolVersion:
    """Container class for outputting version related information from a database."""

    tool_info: str = "x.x.x"
    database_info: str = "x.x.x"


class DatabaseValidationError(Exception):
    """Custom exception for database validation errors."""


class SearchHandler(ABC):
    """
    Abstract class defining tool interface including a database directory / file to search, an input
    query, and an output file to be used for screening.
    """

    def __init__(
        self,
        database_file: str | os.PathLike,
        input_file: str | os.PathLike,
        out_file: str | os.PathLike,
        threads: int = 1,
    ):
        self.db_file = os.path.abspath(database_file)
        self.input_file = os.path.abspath(input_file)
        self.out_file = os.path.abspath(out_file)
        self.threads = threads

        self.arguments_dictionary = {}

        self._validate_db()
        self.version_info = self.get_version_information()

    @property
    def db_directory(self):
        """Directory where databases to be searched are located."""
        return os.path.dirname(self.db_file)

    @property
    def temp_log_file(self):
        """Temporary log file used for this search. Based on outfile name."""
        return f"{self.out_file}.log.tmp"

    @abstractmethod
    def search(self):
        """
        Use a tool to search the input query against a database.
        Should be implemented by all subclasses to perform the actual search against the database.
        """

    @abstractmethod
    def get_version_information(self) -> SearchToolVersion:
        """
        Provide version for the search tool used, to allow reproducibility.
        This method should be implemented by all subclasses to return tool-specific version info.
        """

    def check_output(self):
        """
        Check the output file exists, indicating that the search ran.
        Can be overridden if more complex checks for a particular tool are desired.
        """
        return os.path.isfile(self.out_file)

    def _validate_db(self):
        """
        Validates that the database directory and file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise DatabaseValidationError(
                f"Mandatory screening directory {self.db_directory} not found."
            )

        if not os.path.isfile(self.db_file):
            raise DatabaseValidationError(
                f"Provided database file not found: {self.db_file}."
            )

    @staticmethod
    def is_empty(filepath: str) -> bool:
        """Check if a file is empty or non-existent."""
        try:
            return os.path.getsize(os.path.abspath(os.path.expanduser(filepath))) == 0
        except OSError:
            # Errors such as FileNotFoundError considered empty
            return True

    @staticmethod
    def has_hits(filepath: str) -> bool:
        """Check if a file has any hits (lines that do not start with '#')."""
        try:
            with open(filepath, "r", encoding="utf-8") as file:
                return any(not line.strip().startswith("#") for line in file)
        except FileNotFoundError:
            return False

    def format_args_for_cli(self) -> list:
        """
        Format `self.arguments_dictionary` into a list of strings for use in the command line.
        """
        formatted_args = []
        for key, value in self.arguments_dictionary.items():
            formatted_args.append(str(key))
            if isinstance(value, list):
                formatted_args.append(" ".join(map(str, value)))
            elif value is not None:
                formatted_args.append(str(value))
        return formatted_args

    def run_as_subprocess(self, command, out_file, raise_errors=False):
        """
        Run a command using subprocess.run, piping stdout and stderr to `out_file`.
        """
        logging.debug("SUBPROCESS: %s", " ".join(command))

        with open(out_file, "a", encoding="utf-8") as f:
            result = subprocess.run(
                command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
            )

            if result.returncode != 0:
                command_str = " ".join(command)
                logging.info(
                    "\t ERROR: command %s failed with error %s",
                    command_str,
                    result.stderr,
                )
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                    f" Check {out_file} for logs."
                )

    def __del__(self):
        if os.path.exists(self.temp_log_file):
            os.remove(self.temp_log_file)