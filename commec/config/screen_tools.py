#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
import logging
from commec.config.io_parameters import ScreenIOParameters
from commec.tools.blastn_handler import BlastNHandler
from commec.tools.blastx_handler import BlastXHandler
from commec.tools.blastdmnd_handler import DiamondHandler
from commec.tools.cmscan_handler import CmscanHandler
from commec.tools.hmm_handler import HMMHandler

class ScreenTools():
    """
    Consolidation and initialisation of all databases needed for the commec workflow.
    This could be maybe split into a Biorisk_handler, Taxonomy handler, and Benign handler
    class objects respectively.
    """
    def __init__(self, params : ScreenIOParameters):
        # Biorisk
        self.biorisk_db : HMMHandler = None
        # Taxonomy
        self.protein_db = None
        self.nucleotide_db : BlastNHandler = None
        # Benign
        self.benign_hmm : HMMHandler = None
        self.benign_blastn : BlastNHandler = None
        self.benign_cmscan : CmscanHandler = None

        if params.should_do_biorisk_screening:
            self.biorisk_db = HMMHandler(
                os.path.join(params.db_dir, "biorisk_db"),
                os.path.join(params.db_dir, "biorisk_db/biorisk.hmm"),
                params.query.aa_path,
                f"{params.output_prefix}.biorisk.hmmscan",
            )

        if params.should_do_protein_screening:
            if params.config.search_tool == "blastx":
                self.protein_db = BlastXHandler(
                    os.path.join(params.db_dir, "nr_blast"),
                    os.path.join(params.db_dir, "nr_blast/nr"),
                    input_file = params.query.nt_path,
                    out_file = f"{params.output_prefix}.nr.blastx"
                )
            elif params.config.search_tool == "nr.dmnd" or params.config.search_tool == "diamond":
                self.protein_db = DiamondHandler(
                    os.path.join(params.db_dir, "nr_dmnd"),
                    os.path.join(params.db_dir, "nr_dmnd/nr.dmnd"),
                    input_file = params.query.nt_path,
                    out_file = f"{params.output_prefix}.nr.dmnd"
                )
                if params.config.search_tool == "nr.dmnd":
                    logging.info("""Using old \"nr.dmnd\" keyword for search tool used
                                  will not be supported in future releases, 
                                 consider using \"diamond\" instead.""")
            else:
                raise RuntimeError("Search tool not defined as \"blastx\" or \"diamond\"")

        if params.should_do_nucleotide_screening:
            self.nucleotide_db = BlastNHandler(
                os.path.join(params.db_dir, "nt_blast"),
                os.path.join(params.db_dir, "nt_blast/nt"),
                input_file = f"{params.output_prefix}.noncoding.fasta",
                out_file = f"{params.output_prefix}.nt.blastn"
            )

        if params.should_do_benign_screening:
            self.benign_hmm = HMMHandler(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.hmm"),
                input_file = params.query.nt_path,
                out_file = f"{params.output_prefix}.benign.hmmscan"
            )
            self.benign_blastn = BlastNHandler(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.fasta"),
                input_file = params.query.nt_path,
                out_file = f"{params.output_prefix}.benign.blastn"
            )
            self.benign_cmscan = CmscanHandler(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.cm"),
                input_file = params.query.nt_path,
                out_file = f"{params.output_prefix}.benign.cmscan"
            )
