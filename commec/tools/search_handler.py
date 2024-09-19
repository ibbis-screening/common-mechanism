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
class SearchToolVersion():
    """ Container class for outputting version related information from a database."""
    version_string: str = "x.x.x"
    version_date: str = "Null"
    additional_comment: str = ""

class SearchHandler(ABC):
    """
    Abstract class defining tool interface including a database directory / file to search, an input
    query, and an output file to be used for screening.
    """
    def __init__(self, database_file : str, input_file : str, out_file : str):
        self.db_directory = os.path.dirname(database_file)
        self.db_file = database_file
        self.out_file = out_file
        self.input_file = input_file
        self.temp_log_file = f"{self.out_file}.log.tmp"
        self.arguments_dictionary = None
        self.validate_db()

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

    def validate_db(self):
        """
        Validates that the directory,
        and database file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise FileNotFoundError(f"Database directory not found: {self.db_directory}.")
        if not os.path.isfile(self.db_file):
            raise FileNotFoundError(f"Database file not found: {self.db_file}.")

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

    def get_arguments(self) -> list:
        """
        convert the arguments dictionary into a list,
        structurally ready for appending to a command list of strs.
        """
        my_list = []
        for key, value in self.arguments_dictionary.items():
            my_list.append(str(key))
            if isinstance(value, list):
                my_list.append(" ".join(value))  # Extend the list with all elements in the array
            else:
                my_list.append(str(value))  # Append the value directly if it's not a list
        return my_list

    def run_as_subprocess(self, command, out_file, raise_errors=False):
        """
        Run a command using subprocess.run, piping stdout and stderr to `out_file`.
        """
        logging.debug("SUBPROCESS: %s"," ".join(command))

        with open(out_file, "a", encoding="utf-8") as f:
            result = subprocess.run(
                command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
            )

            if result.returncode != 0:
                command_str = ' '.join(command)
                logging.info("\t ERROR: command %s failed with error %s", command_str, result.stderr)
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                    f" Check {out_file} for logs."
                )
    def __del__(self):
        if os.path.exists(self.temp_log_file):
            os.remove(self.temp_log_file)