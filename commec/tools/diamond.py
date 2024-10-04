#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Handler for Diamond blastx search of protein databases using nucleotide query.
Instantiate a DiamondHnadler, with input local database, input fasta, and output file.
Alter member variables after instantiation to change behaviour, before calling search().
Throws if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""
import os
import glob
import subprocess
from typing import Optional
import logging
from multiprocessing import Pool

from commec.tools.blast_tools import BlastHandler
from commec.tools.search_handler import SearchToolVersion


class DiamondHandler(BlastHandler):
    """
    Diamond blastx search handler. Initialise instance with inputs and outputs,
    alter arguments as desired, and call search() to run Diamond blastx via CLI.
    Concatenates all diamond outputs into a single output file.
    """

    def __init__(self, database_file: str, input_file: str, out_file: str):
        super().__init__(database_file, input_file, out_file)
        self.frameshift: int = 15
        self.do_range_culling = True
        self.threads = 1
        self.jobs: Optional[int] = None
        self.output_format = "6"
        self.output_format_tokens = [
            "qseqid",
            "stitle",
            "sseqid",
            "staxids",
            "evalue",
            "bitscore",
            "pident",
            "qlen",
            "qstart",
            "qend",
            "slen",
            "sstart",
            "send",
        ]

    def run_as_subprocess_args(self, args):
        """Wrapper for run_as_subprocess, for Diamond pooling processes."""
        command, out_file = args
        self.run_as_subprocess(command, out_file, True)

    def determine_runs_threads_and_cycles(
        self, max_threads: int, n_database_files: int
    ):
        """
        Determine how many diamonds runs, and how many processors per run, need to occur
        for best CPU use and efficiency.
        """
        number_of_databases: int = n_database_files
        highest_common_denominator: int = max(number_of_databases, 1)
        while highest_common_denominator > 1:
            if (
                number_of_databases % highest_common_denominator == 0
                and max_threads % highest_common_denominator == 0
            ):
                break
            highest_common_denominator -= 1

        n_concurrent_diamond_runs: int = (
            highest_common_denominator if self.jobs is None else self.jobs
        )
        n_threads_per_diamond_run: int = int(max_threads) / int(
            n_concurrent_diamond_runs
        )
        n_cycles: int = int(number_of_databases) / int(n_concurrent_diamond_runs)
        return n_concurrent_diamond_runs, n_threads_per_diamond_run, n_cycles

    def search(self):
        # Find all files matching the pattern nr*.dmnd in DB_PATH
        db_files = glob.glob(f"{self.db_directory}/nr*.dmnd")
        max_threads: int = self.threads

        # Sanity checks on job and thread settings.
        if len(db_files) == 0:
            raise FileNotFoundError(
                f"Mandatory Diamond database directory {self.db_directory} "
                "contains no databases!"
            )

        print(self.threads)
        print(len(db_files))
        print(self.jobs)
        n_concurrent_diamond_runs, n_threads_per_diamond_run, n_cycles = (
            self.determine_runs_threads_and_cycles(self.threads, len(db_files))
        )

        logging.info(
            "Using %i concurrent diamond runs, with %i processes per run, across %i cycles.",
            n_concurrent_diamond_runs,
            n_threads_per_diamond_run,
            n_cycles,
        )

        # We keep the below checks, incase the user inputs a -j argument, that is a poor choice.
        if n_concurrent_diamond_runs < 1:
            n_concurrent_diamond_runs = 1
            logging.info(
                "WARNING, numnber of Diamond job pools cannot be < 1. "
                "Number of job pools as been reset to 1."
            )
        if n_threads_per_diamond_run < 1:
            n_threads_per_diamond_run = 1
            logging.info(
                "WARNING, number of indiviudal Diamond job allocated threads cannot be < 1. "
                "The number of threads per Diamond job as been reset to 1."
            )
        if n_threads_per_diamond_run * n_concurrent_diamond_runs < max_threads:
            logging.info(
                "WARNING, total number of threads across concurrent Diamond job pools [%i*%i]"
                " is less than number of allocated threads[%i]. CPU may be underutilised.",
                n_threads_per_diamond_run,
                n_concurrent_diamond_runs,
                max_threads,
            )

        if n_threads_per_diamond_run * n_concurrent_diamond_runs > max_threads:
            logging.info(
                "WARNING, number of Diamond job pools [%i] are using [%i] threads"
                " each. CPU may be bottlenecked.",
                n_concurrent_diamond_runs,
                n_threads_per_diamond_run,
            )
        if len(db_files) < n_concurrent_diamond_runs:
            logging.info(
                "WARNING, The Diamond database has only been split into %i parts."
                " Excessive number of requested concurrent Diamond jobs %i",
                len(db_files),
                n_concurrent_diamond_runs,
            )
        if len(db_files) % n_concurrent_diamond_runs > 0:
            logging.info(
                "WARNING, The number of Diamond database files [%i] is not "
                " divisible by concurrent jobs [%i]. CPU may be underutilised.",
                len(db_files),
                n_concurrent_diamond_runs,
            )

        # Lists to track each Diamond: command inputs, outputs, and file logs.
        pool_arguments = []
        output_files = []
        output_lognames = []

        # Generate an Array of tuple inputs for passing as arguments to Pool.
        for i, db_file in enumerate(db_files, 1):
            # Generate Diamond command string:
            output_file = f"{self.out_file}.{i}.tsv"
            output_files.append(output_file)
            command = [
                "diamond",
                "blastx",
                "--quiet",
                "-d",
                db_file,
                "--threads",
                str(n_threads_per_diamond_run),
                "-q",
                self.input_file,
                "-o",
                output_file,
                "--frameshift",
                str(self.frameshift),
                "--outfmt",
                self.output_format,
            ]
            command.extend(self.output_format_tokens)
            if self.do_range_culling:
                command.append("--range-culling")

            # Generate log filenames:
            log_i_filename: str = self.temp_log_file + "_" + str(i)
            output_lognames.append(log_i_filename)

            pool_arguments.append((command, log_i_filename))

        # Run each input, with pooling.
        with Pool(n_concurrent_diamond_runs) as pool:
            pool.map(self.run_as_subprocess_args, pool_arguments)

        # Concatenate all output files, and remove the temporary ones.
        with open(self.out_file, "w", encoding="utf-8") as outfile:
            for o_file in output_files:
                if os.path.exists(o_file):
                    with open(o_file, "r", encoding="utf-8") as infile:
                        outfile.write(infile.read())
                    os.remove(o_file)

        # Logs are concatenated into a global diamond log.
        with open(self.temp_log_file, "w", encoding="utf-8") as outlog:
            for log_f in output_lognames:
                if os.path.exists(log_f):
                    with open(log_f, "r", encoding="utf-8") as infile:
                        outlog.write(infile.read())
                    os.remove(log_f)

    def get_version_information(self) -> SearchToolVersion:
        try:
            db_files = glob.glob(f"{self.db_file}*dmnd")
            result = subprocess.run(
                ["diamond", "dbinfo", "--db", db_files[0]],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.splitlines()
            database_info: str = lines[1]

            result = subprocess.run(
                ["diamond", "version"], capture_output=True, text=True, check=True
            )

            tool_info: str = result.stdout.strip()
            return SearchToolVersion(tool_info, database_info)
        except subprocess.CalledProcessError:
            return None
