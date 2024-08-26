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

from commec.io_parameters import ScreenIOParameters

class CommecDatabases():
    """
    Consolidation and initialisation of all databases needed for the commec workflow.
    This could be maybe split into a Biorisk_handler, Taxonomy handler, and Benign handler
    class objects respectively.
    """
    def __init__(self, params : ScreenIOParameters):
        # Biorisk
        self.biorisk_db : HMMDataBase = None
        # Taxonomy
        self.protein_db : BlastXDataBase = None
        self.nucleotide_db : BlastNDataBase = None
        # Benign
        self.benign_hmm : HMMDataBase = None
        self.benign_blastn : BlastNDataBase = None
        self.benign_cmscan : CmscanDataBase = None

        if params.should_do_biorisk_screening:
            self.biorisk_db = HMMDataBase(
                os.path.join(params.db_dir, "biorisk_db"),
                os.path.join(params.db_dir, "biorisk_db/biorisk.hmm"),
                params.query.fasta_aa_filepath,
                f"{params.output_prefix}.biorisk.hmmscan",
            )

        if params.should_do_protein_screening:
            if params.inputs.search_tool == "blastx":
                self.protein_db = BlastXDataBase(
                    os.path.join(params.db_dir, "nr_blast"),
                    os.path.join(params.db_dir, "nr_blast/nr"),
                    input_file = params.query.fasta_aa_filepath,
                    out_file = f"{params.output_prefix}.nr.blastx"
                )
            elif params.inputs.search_tool == "nr.dmnd":
                self.protein_db = BlastXDataBase(
                    os.path.join(params.db_dir, "nr_dmnd"),
                    os.path.join(params.db_dir, "nr_blast/nr.dmnd"),
                    input_file = params.query.fasta_aa_filepath,
                    out_file = f"{params.output_prefix}.nr.dmnd"
                )
            else:
                raise RuntimeError("Search tool not defined as \"blastx\" or \"nr.dmnd\"")

        if params.should_do_nucleotide_screening:
            self.nucleotide_db = BlastNDataBase(
                os.path.join(params.db_dir, "nt_blast"),
                os.path.join(params.db_dir, "nt_blast/nt"),
                input_file = f"{params.output_prefix}.noncoding.fasta",
                out_file = f"{params.output_prefix}.nt.blastn"
            )

        if params.should_do_benign_screening:
            self.benign_hmm = HMMDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.hmm"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.hmmscan"
            )
            self.benign_blastn = BlastNDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.fasta"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.blastn"
            )
            self.benign_cmscan = CmscanDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.cm"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.cmscan"
            )

@dataclass
class DatabaseVersion():
    """ Container class for outputting version related information from a database."""
    def __init__(self, input_version : str = "x.x.x", input_date : str = "Null", input_comment :str = ""):
        self.version_string = input_version
        self.version_date = input_date
        self.additional_comment = input_comment
    version_string : str
    version_date : str
    additional_comment : str

class DatabaseHandler():
    """ 
    Abstract Class for holding the directory, and file of a database, 
    including input and output files for screening, as well as screening arguments. 
    Override the Screen() method for a custom database implementation.
    """
    # Start Database Handler API
    def screen(self):
        """ Virtual function to be called by any child database implementation"""
        raise NotImplementedError(
            "This class must override the .screen() to use DatabaseHandler."
        )

    def get_version_information(self) -> DatabaseVersion:
        """ Override to call version retrieval information on a database."""
        new_version_info = DatabaseVersion("", "", "Version retrieval not yet supported for this database")
        return new_version_info

    def check_output(self):
        """ 
        Simply checks for the existance of the output file, 
        indicating that the database screening ran. Can be overrided if
        more complex checks for a particular db is desired.
        """
        if os.path.isfile(self.out_file): # Consider adding checks from File_Tools such that empty outputs are invalid?
            return True
        return False

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        self.db_directory = directory
        self.db_file = database_file
        self.out_file = out_file
        self.input_file = input_file
        self.temp_log_file = f"{self.out_file}.log.tmp"
        self.arguments_dictionary = {}
        self.validate_directory()

    # End override API

    def check_input(self):
        """ 
        Simply checks for the existance of the input file, This is separated 
        from database location, as sometimes the input file is generated at a 
        previous step, and may not exist during DatabaseHandler instantiation.
        """
        if os.path.isfile(self.input_file):
            return True
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

    def validate_directory(self):
        """ 
        Validates that the directory, 
        and database file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise FileNotFoundError(f"Mandatory screening directory {self.db_directory} not found.")
        if not os.path.isfile(self.db_file):
            logging.error("Bad database file!")
            #raise FileNotFoundError(f"Mandatory screening directory {self.db_file} not found.")

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
            return
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

    def get_version_information(self) -> DatabaseVersion:
        """ At the moment this is just grabbing some basic header info of the 
        first entrant of the hmm database. Not really a true version control.
        But better than nothing at the moment. There may be some way to return
        some version information from hmmcan itself, look into that. """
        version : str = None
        date : str = None
        comment : str = None
        try:
            with open(self.db_file, 'r', encoding = "utf-8") as file:
                for line in file:
                    if line.startswith("NAME"):
                        comment = line.split(maxsplit=1)[1].strip()
                        continue
                    if line.startswith("DATE"):
                        date = line.split(maxsplit=1)[1].strip()
                        continue
                    if line.startswith("HMMER3/f"):
                        version = line.split("[",maxsplit=1)[1].split("|")[0].strip()
                        continue
                    # Early exit if data has been found
                    if version and date and comment:
                        break
            return DatabaseVersion(version, date, comment)

        except RuntimeError:
            return super().get_version_information()

class CmscanDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Hmmer files for commec screening. """
    def screen(self):
        command = [
            "cmscan", 
            "--tblout",
            self.out_file,
            self.db_file,
            self.input_file
            ]
        self.run_as_subprocess(command, self.temp_log_file)

class BlastXDataBase(DatabaseHandler):
    """ A Database handler specifically for use with BlastX files for commec screening. 
    Allows for full customization of any of the callable blast flags. A better implementation
    may be achieved using import blast"""

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
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
    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
        # We fill this with the defaults, however they can always be overridden after instancing, before screening.
        self.arguments_dictionary = {
            "-num_threads": 8,
            "-evalue": "10",
            "-max_target_seqs": 50,
            "-culling_limit": 5,
            "-outfmt": ["7", "qacc",   "stitle",   "sacc",   "staxids",
                        "evalue", "bitscore", "pident", "qlen",
                        "qstart", "qend",     "slen",   "sstart",  "send"]
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

#TODO: Ensure updates from Tessas handling of Diamond are correctly incorperated here.
class DiamondDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Diamond files for commec screening. """

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        super().__init__(directory, database_file, input_file, out_file)
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