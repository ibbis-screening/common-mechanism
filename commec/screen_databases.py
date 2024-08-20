#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
import glob
from dataclasses import dataclass
import subprocess
import logging

@dataclass
class DatabaseHandler():
    """ 
    Abstract Class for holding the directory, and file of a database, 
    including input and output files for screening, as well as screening arguments. 
    Override the Screen() method for a custom database implementation.
    """
    # Start Database Handler API
    def screen(self):
        """ Virtual function to be called by any child database implementation"""
        raise NotImplementedError("This class must override the .screen() to use DatabaseHandler.")

    def check_output(self):
        """ 
        Simply checks for the existance of the output file, 
        indicating that the database screening ran. Can be overrided if
        more complex checks for a particular db is desired.
        """
        if os.path.isfile(self.out_file): # Consider adding checks from File_Tools such that empty outputs are invalid?
            return True
        return False

    def __init__(self, directory : str, databse_file : str, input_file : str, out_file : str):
        self.db_directory = directory
        self.db_file = databse_file
        self.out_file = out_file
        self.input_file = input_file
        self.temp_log_file = f"{self.input_file}.log.tmp"
        self.arguments_dictionary = {}
        self.validate_directory()

    # End override API
    def get_arguments(self) -> list:
        """ 
        convert the arguments dictionary into a list, 
        structurally ready for appending to a command list of strs.
        """
        my_list = []
        for key, value in self.arguments_dictionary.items():
            my_list.append(str(key))
            if isinstance(value, list):
                my_list.extend(str(v) for v in value)  # Extend the list with all elements in the array
            else:
                my_list.append(str(value))  # Append the value directly if it's not a list
        return my_list

    def validate_directory(self):
        """ 
        Validates that the directory, 
        and database file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise FileNotFoundError(f"Mandatory screening directory {self.db_directory} not found.")
        if not os.path.isfile(self.db_file):
            raise FileNotFoundError(f"Mandatory screening directory {self.db_file} not found.")

    def run_as_subprocess(self, command, out_file, raise_errors=False):
        """
        Run a command using subprocess.run, piping stdout and stderr to `out_file`.
        """
        with open(out_file, "a", encoding="utf-8") as f:
            result = subprocess.run(
                command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
            )
            if result.returncode != 0:
                command_str = ' '.join(command)
                logging.info("\t ERROR: command %s failed", command_str)
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                    f" Check {out_file} for logs."
                )
            
    def __del__(self):
        if os.path.exists(self.temp_log_file):
            os.remove(self.temp_log_file)

class HMMDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Hmmer files for commec screening. """
    def screen(self):
        command = [
            "hmmscan", 
            "--domtblout",
            self.out_file,
            self.db_file,
            self.input_file
            ]
        self.run_as_subprocess(command, self.temp_log_file)

class BlastXDataBase(DatabaseHandler):
    """ A Database handler specifically for use with BlastX files for commec screening. 
    Allows for full customization of any of the callable blast flags. A better implementation
    may be achieved using import blast"""

    def __init__(self, directory : str, databse_file : str, input_file : str, out_file : str):
        super().__init__(directory, databse_file, input_file, out_file)
        # We fill this with the defaults, however they can always be overridden after instancing, before screening.
        self.arguments_dictionary = {
            "-num_threads": 1,
            "-evalue": "1e-10",
            "-word_size": 6,
            "-threshold": 21,
            "-max_target_seqs": 5000,
            "-culling_limit": 50,
            "-window_size": 40,
            "-matrix": "BLOSUM62",
            "-gapopen": 11,
            "-gapextend": 1,
            "-seg": "yes",
            "-outfmt": ["7", "qacc",   "stitle",   "sacc",   "staxids",
                        "evalue", "bitscore", "pident", "qlen",
                        "qstart", "qend",     "slen",   "sstart",  "send"]
        }
        self.blastcall = "blastx"

    def screen(self):
        command = [
            self.blastcall,
            "-db", self.db_file,
            "-query", self.input_file,
            "-out", self.out_file,
        ]
        command.extend(self.get_arguments())
        self.run_as_subprocess(command,self.temp_log_file)

class BlastNDataBase(DatabaseHandler):
    """ A Database handler specifically for use with BlastX files for commec screening. 
    Allows for full customization of any of the callable blast flags. A better implementation
    may be achieved using import blast"""
    def __init__(self, directory : str, databse_file : str, input_file : str, out_file : str):
        super().__init__(directory, databse_file, input_file, out_file)
        # We fill this with the defaults, however they can always be overridden after instancing, before screening.
        self.arguments_dictionary = {
            "-num_threads": 8,
            "-evalue": "10",
            "-max_target_seqs": 50,
            "-culling_limit": 5,
            "-outfmt": ["7", "qacc", "stitle", "sacc", "staxids",
                        "evalue", "bitscore", "pident", "qlen",
                        "qstart", "qend", "slen", "sstart", "send"]
        }
        self.blastcall = "blastn"

    def screen(self):
        command = [
            self.blastcall,
            "-db", self.db_file,
            "-query", self.input_file,
            "-out", self.out_file,
        ]
        command.extend(self.get_arguments())
        self.run_as_subprocess(command,self.temp_log_file)

#TODO: Await updates from Tessas handling of Diamond, to incorperate them here.
class DiamondDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Diamond files for commec screening. """

    def __init__(self, directory : str, databse_file : str, input_file : str, out_file : str):
        super().__init__(directory, databse_file, input_file, out_file)
        self.frameshift = 15
        self.do_range_culling = True
        self.threads = 1
        self.jobs = None
        self.output_format = "6"
        self.output_format_tokens = [
            "qseqid",   "stitle",   "sseqid",   "staxids",
            "evalue", "bitscore", "pident", "qlen",
            "qstart", "qend",     "slen",   "sstart",  "send"]

    def get_output_format(self) -> str:
        """ Returns a formatted string version of the output format for blastx"""
        return self.output_format + " ".join(self.output_format_tokens)

    def screen(self):
        # Find all files matching the pattern nr*.dmnd in DB_PATH
        db_files = glob.glob(f"{self.db_directory}/nr*.dmnd")

        if len(db_files) == 0:
            raise FileNotFoundError(f"Mandatory Diamond database directory {self.db_directory} contains no databases!")

        # Run diamond blastx in parallel using subprocess
        processes = []
        for i, db_file in enumerate(db_files, 1):
            output_file = f"{self.out_file}.{i}.tsv"
            command = [
                "diamond", "blastx",
                "--quiet",
                "-d", db_file,
                "--threads", str(self.threads),
                "-q", self.input_file,
                "-o", output_file,
                "--outfmt", self.get_output_format(),
                "--frameshift", str(self.frameshift),
            ]

            if self.do_range_culling:
                command.append("--range-culling")

            if self.jobs is not None:
                command.extend(["-j", str(self.jobs)])

            process = subprocess.Popen(command, self.temp_log_file)
            processes.append(process)

        # Wait for all processes to finish
        for process in processes:
            process.wait()

        # Concatenate all output files
        with open("self.out_file", 'w', encoding="utf-8") as outfile:
            for i in range(1, len(db_files) + 1):
                with open(f"{self.out_file}.{i}.tsv", 'r', encoding="utf-8") as infile:
                    outfile.write(infile.read())

        # Remove the individual output files
        for i in range(1, len(db_files) + 1):
            os.remove(f"{self.out_file}.{i}.tsv")