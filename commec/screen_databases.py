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
        self.tax_dir = os.path.join(database_dir, "taxonomy")
        self.nt_dir = os.path.join(database_dir, "nt_blast")

        self.biorisk_db = os.path.join(self.biorisk_dir, "biorisk.hmm")
        self.benign_db = os.path.join(self.benign_dir, "benign.hmm")
        self.nt_db = os.path.join(self.nt_dir, "nt")


    @property
    def nr_dir(self):
        """
        Get the directory containing the needed NR database (which depends on the tool used).
        """
        if self.search_tool == "blastx":
            return os.path.join(self.dir, "nr_blast")
        elif self.search_tool == "diamond":
            return os.path.join(self.dir, "nr_dmnd")
        else:
            raise ValueError(
                "Protein search tool must be either 'blastx' or 'diamond'."
            )
    
    @property
    def nr_db(self):
        if self.search_tool == "blastx":
            return os.path.join(self.nr_dir, "swissprot")
        elif self.search_tool == "diamond":
            return os.path.join(self.nr_dir, "nr.dmnd")
        else:
            raise ValueError(
                "Protein search tool must be either 'blastx' or 'diamond'."
            )

    def validate(self, in_fast_mode=False, skip_nt_search=False):
        """
        Make sure all the needed databases exist.
        """
        needed_dirs = [self.biorisk_dir, self.benign_dir] 
        # Fast mode doesn't need to check protein search directories
        if not in_fast_mode:
            needed_dirs.append(self.nr_dir)
            needed_dirs.append(self.tax_dir)
        
        if not in_fast_mode and not skip_nt_search:
            needed_dirs.append(self.nt_dir)

        for db_dir in needed_dirs:
            if not os.path.isdir(db_dir):
                raise FileNotFoundError(
                    f"Mandatory screening directory {db_dir} not found."
                )

        for db in [self.biorisk_db, self.benign_db]:
            if not os.path.isfile(db):
                raise FileNotFoundError(f"Mandatory screening database {db} not found.")
