#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `ScreenDatabases` class.
"""
import os

class ScreenDatabases:
    """
    The Common Mechanism depends on a lot of datbases in a particular structure (encoded here).
    """

    def __init__(self, database_dir, search_tool):
        self.dir = database_dir
        self.search_tool = search_tool
        self.biorisk_dir = os.path.join(database_dir, "biorisk_db")
        self.benign_dir = os.path.join(database_dir, "benign_db")
        self.nt_dir = os.path.join(database_dir, "nt_blast/nt")

        self.biorisk_db = os.path.join(self.biorisk_dir, "biorisk.hmm")
        self.benign_db = os.path.join(self.benign_dir, "benign.hmm")

    @property
    def nr_dir(self):
        """
        Get the directory containing the needed NR database (which depends on the tool used).
        """
        if self.search_tool == "blastx":
            return os.path.join(self.dir, "nr_blast/nr")
        elif self.search_tool == "diamond":
            return os.path.join(self.dir, "nr_dmnd")
        else:
            raise ValueError(
                "Protein search tool must be either 'blastx' or 'diamond'."
            )

    def validate(self, in_fast_mode=False):
        """
        Make sure all the needed databases exist.
        """
        for db_dir in [self.biorisk_dir, self.benign_dir]:
            if not os.path.isdir(db_dir):
                raise FileNotFoundError(
                    f"Mandatory screening directory {db_dir} not found."
                )

        for db in [self.biorisk_db, self.benign_db]:
            if not os.path.isfile(db):
                raise FileNotFoundError(f"Mandatory screening database {db} not found.")

        # Fast mode doesn't need to check protein search directories
        if not in_fast_mode and not os.path.isdir(self.nr_dir):
            raise FileNotFoundError(
                f"Protein screening directory {self.nr_dir} not found."
            )
