#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science 
"""
Module for CLI setup of Commec, such that required 
databases are downloaded in a desired database directory.
"""
import sys
import os
import argparse
import subprocess
import urllib.request
import zipfile

DESCRIPTION = """Helper script for downloading the databases
 required for running the Common Mechanism Screen"""

class CliSetup:
    """
    Interacts with the user via the CLI, and then downloads appropriate files
    into appropraite folders as requested by the user, to allow for the Commec
    Screen workflow.

    The current decision flow is the following:
        > Start
            > Specify overall database directory
            > Confirm Biorisk download URL
            > decide Blast protein database
                > Choose Blast protein database
            > decide Blast Nucleotide database
                > Choose Blast Nucleotide database
            > decide Taxonomy database
                > Confirm taxonomy URL
            > Confirm.
        
        Going "back" takes you back up the tree.
        Deciding not to something skips the trees children.

        Due to a lack of wanting to over-engineer this, 
        each step tree is implemented using a function
        collectively acting as a pseudo-statemachine.
        Once Start is successfully returned, the setup process can begin.
    """

    def __init__(self, automate : bool = False):
        self.ask : bool = not automate
        self.database_directory : str = "./commec-dbs/"

        self.download_biorisk : bool = True
        self.default_biorisk_download_url : str = "https://f005.backblazeb2.com/file/common-mechanism-dbs/common-mechanism-dbs.zip"
        self.biorisk_download_url : str = self.default_biorisk_download_url
        
        self.download_blastnr : bool = False
        self.blastnr_database : str = "nr"
        self.download_blastnt : bool = False
        self.blastnt_database : str = "nt"

        self.download_example_blastnr : bool = False
        self.download_example_blastnt : bool = False

        self.download_taxonomy : bool = False
        self.default_taxonomy_download_url : str = "ftp\://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
        self.taxonomy_download_url : str = self.default_taxonomy_download_url

        # Hard coded, sure, but gives a list of valid dbs to check against.
        # A smarter implementation would take the time to call the --showall function without pretty,
        # and then populate this array with the outputs.
        # At least this way we can slightly curate out the non-aa and non-nt ones.
        self.protein_database_names = [
            "BLASTDB", "Betacoronavirus", "28S_fungal_sequences", "16S_ribosomal_RNA", 
            "18S_fungal_sequences", "ITS_RefSeq_Fungi", "ITS_eukaryote_sequences", 
            "LSU_eukaryote_rRNA", "LSU_prokaryote_rRNA", "SSU_eukaryote_rRNA", 
            "env_nt", "env_nr", "human_genome", "landmark", "mito", "mouse_genome", 
            "nr", "nt_euk", "nt", "nt_others", "nt_prok", "nt_viruses", "pataa", 
            "patnt", "pdbaa", "pdbnt", "ref_euk_rep_genomes", "ref_prok_rep_genomes", 
            "ref_viroids_rep_genomes", "ref_viruses_rep_genomes", "refseq_select_rna", 
            "refseq_select_prot", "refseq_protein", "refseq_rna", "swissprot", "tsa_nr", 
            "tsa_nt", "taxdb", "core_nt"
        ]

        self.nucleotide_database_names = [
            "BLASTDB", "Betacoronavirus", "28S_fungal_sequences", "16S_ribosomal_RNA", 
            "18S_fungal_sequences", "ITS_RefSeq_Fungi", "ITS_eukaryote_sequences", 
            "LSU_eukaryote_rRNA", "LSU_prokaryote_rRNA", "SSU_eukaryote_rRNA", 
            "env_nt", "env_nr", "human_genome", "landmark", "mito", "mouse_genome", 
            "nr", "nt_euk", "nt", "nt_others", "nt_prok", "nt_viruses", "pataa", 
            "patnt", "pdbaa", "pdbnt", "ref_euk_rep_genomes", "ref_prok_rep_genomes", 
            "ref_viroids_rep_genomes", "ref_viruses_rep_genomes", "refseq_select_rna", 
            "refseq_select_prot", "refseq_protein", "refseq_rna", "swissprot", "tsa_nr", 
            "tsa_nt", "taxdb", "core_nt"
        ]

        if automate:
            self.do_setup()
        else:
            self.start()

    def start(self):
        """
        Starts the user interrogation process.
        """
        print("""                       Welcome to
 ██████╗ ██████╗ ███╗   ███╗███╗   ███╗███████╗ ██████╗
██╔════╝██╔═══██╗████╗ ████║████╗ ████║██╔════╝██╔════╝
██║     ██║   ██║██╔████╔██║██╔████╔██║█████╗  ██║     
██║     ██║   ██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██║     
╚██████╗╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗╚██████╗
 ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝ ╚═════╝""")
        print("                 The Common Mechanism!",
              "\n\nCopyright © 2021-2024 International Biosecurity",
              "\nand Biosafety Initiative for Science",
              "\n\nThis script will help download the mandatory databases ",
              "\nrequired for using Commec Screen, and requires a stable",
              "\ninternet connection, wget, and update_blastdb.pl.",
              "\n -> You can exit this setup at any time with \"exit\"",
              "\n -> You can return to a previous step with \"back\"\n")
        
        self.setup_overall_directory()
        self.do_setup()
        
    def setup_overall_directory(self):
        """
        Get user inputs for global directory to store databases.
        """
        self.print_step(1)
        user_input : str = ""
        print("\nPlease provide the absolute or relative filepath",
              "to where you would like the Commec databases to be located...",
              "\nPress <Enter> to use default: ",
              self.database_directory)
        while(True):
            user_input : str = self.user_input()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.stop() # This is the first step, no point going back.
            if len(user_input) > 0:
                self.database_directory = user_input

            # Consider adding some sort of check here, that the user
            # input is mkdir compatible.
            try:
                subprocess.run(["mkdir","-p", self.database_directory], check = True)
                subprocess.run(["rm","-rf",self.database_directory], check = True)
            except subprocess.CalledProcessError:
                print(user_input, " is not a valid directory structure!")
                continue

            print("Using database directory: ", self.database_directory)
            self.get_biorisk_url()
            return

    def get_biorisk_url(self):
        """
        Get the URL where the Commec Biorisk and Benign databases are located.
        """
        self.print_step(2)
        user_input : str = ""
        print("Please provide the URL to download the Commec database.",
              "\nPress <Enter> to use existing: ",
              self.biorisk_download_url)
        while(True):
            user_input : str = self.user_input()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.setup_overall_directory()
                return
            if len(user_input) > 0:
                self.biorisk_download_url = user_input

            # Consider adding some sort of check here, that the user
            # input is mkdir compatible.
            print("Checking URL is valid ... ")
            if not self.check_url_exists(self.biorisk_download_url):
                print(self.biorisk_download_url, " is not a valid URL! (or you are not connected to the internet)")
                continue
            
            print("Using Commec Biorisk and Benign URL:", self.biorisk_download_url)
            self.decide_blastnr()
            return
    
    def decide_blastnr(self):
        """ Decide whether a Nucleotide database needs to be downloaded. """
        self.print_step(3,1)
        print("Do you want to download a databases for protein screening?",
              "\n\"y\" or \"n\", for yes or no.")
        while(True):
            user_input : str = input(">>>").strip().lower()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.get_biorisk_url()
                return
            if user_input == "y" or user_input == "yes":
                self.download_blastnr = True
                self.get_blastnr()
                return
            if user_input == "n" or user_input == "no":
                self.download_blastnr = False
                self.decide_blastnt()
                return
            print("Unrecognised input (", user_input, ")")

    def get_blastnr(self):
        """ Decide what Protein database needs to be downloaded. """
        self.print_step(3,2)
        print("Which database should be used for taxonomy protein screening?",
              "\ntype \"options\" to fetch the list of options.",
              "\nPress <Enter> to use existing: ",
              self.blastnr_database)
        while(True):
            user_input : str = input(">>>").strip().lower()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.decide_blastnr()
                return
            if user_input == "options":
                self.print_database_options()
                continue
            if len(user_input) == 0:
                print("Using NR database:", self.blastnr_database)
                self.decide_blastnt()
                return
            if user_input in self.protein_database_names:
                self.blastnr_database = user_input
                print("Using NR database:", self.blastnr_database)
                self.decide_blastnt()
                return
            print("Unrecognised or invalid input (", user_input, ")")

    def decide_blastnt(self):
        """ Decide what Protein database needs to be downloaded. """
        self.print_step(4,1)
        print("Do you want to download a databases for nucleotide screening?",
              "\n\"y\" or \"n\", for yes or no.")
        while(True):
            user_input : str = input(">>>").strip().lower()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.get_biorisk_url()
                return
            if user_input == "y" or user_input == "yes":
                self.download_blastnt = True
                self.get_blastnt()
                return
            if user_input == "n" or user_input == "no":
                self.download_blastnt = False
                self.decide_taxonomy()
                return
            print("Unrecognised input (", user_input, ")")

    def get_blastnt(self):
        """ Decide what Nucleotide database needs to be downloaded. """
        self.print_step(4,2)
        print("Which database should be used for taxonomy nucleotide screening?",
              "\ntype \"options\" to fetch the list of options.",
              "\nPress <Enter> to use existing: ",
              self.blastnt_database)
        while(True):
            user_input : str = input(">>>").strip().lower()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.decide_blastnr()
                return
            if user_input == "options":
                self.print_database_options()
                continue
            if len(user_input) == 0:
                print("Using NT database:", self.blastnt_database)
                self.decide_taxonomy()
                return
            if user_input in self.nucleotide_database_names:
                self.blastnt_database = user_input
                print("Using NT database:", self.blastnt_database)
                self.decide_taxonomy()
                return
            print("Unrecognised or invalid input (", user_input, ")")

    def decide_taxonomy(self):
        """ Decide whether taxonomy database need to be downloaded. """
        self.print_step(5,1)
        print("Do you want to download the Taxonomy databases? ( less than 1 GB)",
              "\n\"y\" or \"n\", for yes or no.")
        while(True):
            user_input : str = input(">>>").strip().lower()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.decide_blastnt()
                return
            if user_input == "y" or user_input == "yes":
                self.download_taxonomy = True
                self.confirm()
                return
            if user_input == "n" or user_input == "no":
                self.download_taxonomy = False
                self.confirm()
                return
            print("Unrecognised input (", user_input, ")")

    def get_taxonomy_url(self):
        """
        DEPRECATED - Easier to get this with less bugs through update_blastdbs.pl
        Get the URL where the Taxonomy database is located.
        """
        self.print_step(5,2)
        user_input : str = ""
        print("Please provide the URL to download the BIORISK database.",
              "\nPress <Enter>, to use existing: ",
              self.taxonomy_download_url)
        while(True):
            user_input : str = self.user_input()
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.setup_overall_directory()
                return
            if len(user_input) > 0:
                self.taxonomy_download_url = user_input

            # Consider adding some sort of check here, that the user
            # input is mkdir compatible.
            print("Checking URL is valid ... ")
            if not self.check_url_exists(self.taxonomy_download_url):
                print(self.taxonomy_download_url, " is not a valid URL! (or you are not connected to the internet)")
                continue
            
            print("Using Biorisk URL: %s", self.database_directory)
            self.confirm()
            return
        
    def confirm(self):
        """ Simply allows the user one last chance to confirm their settings."""
        self.print_step(6)
        print("The following settings will be used to setup Commec:",
              "\n -> Database Directory: ", self.database_directory,
              "\n -> Commec Biorisk and Benign URL: ", self.biorisk_download_url)
        if (self.download_blastnr):
            print(" -> PROTEIN NR Database: ", self.blastnr_database)
        if (self.download_blastnt):
            print(" -> NUCLEOTIDE NT Database: ", self.blastnr_database)
        if(self.download_taxonomy):
            print(" -> Taxonomy database will be downloaded.")

        print("\n\nPress <Enter> to confirm these settings, ",
              "\ntype \"back\" to alter previous settings, ",
              "\ntype \"exit\" to abort setup.",
              "\ntype \"restart\" to go back to the beginning.")
        while(True):
            user_input : str = input("").strip().lower()
            if len(user_input) == 0:
                return
            if user_input == "exit":
                self.stop()
            if user_input == "back":
                self.get_biorisk_url()
                return
            if user_input == "restart":
                self.setup_overall_directory()
                return
            print("Unrecognised input (", user_input, ")")

    def do_setup(self):
        """ 
        Once CliSetup state has been finalized,
        call to perform the required actions. 
        """
        self.print_step(7)

        os.makedirs(self.database_directory, exist_ok=True)

        if self.download_biorisk:
            command = ["wget","-c","-P",self.database_directory, self.biorisk_download_url]
            print(
                "Downloading Biorisk database from\n", 
                self.biorisk_download_url
            )
            result = subprocess.run(command, check = True)
            if result.returncode != 0:
                command_str = " ".join(command)
                print(
                    "\t ERROR: Command",
                    command_str,
                    "failed with error",
                    result.stderr,
                )
            # Parse the URL to extract the path
            parsed_url = urllib.parse.urlparse(self.biorisk_download_url)
            filename_zipped = os.path.join(self.database_directory,os.path.basename(parsed_url.path))

            print("Extracting Biorisk databases...")
            # Open the zip file and extract its contents
            with zipfile.ZipFile(filename_zipped, 'r') as zip_ref:
                zip_ref.extractall(self.database_directory)

            #subprocess.run(["rm","-rf", filename_zipped], check = True)
            os.remove(filename_zipped)

        if self.download_blastnr:
            print(
                "Downloading ",
                self.blastnr_database,
                "Blast database for Protein Screening from NCBI \n"
            )
            nr_directory = os.path.join(self.database_directory,"nr_blast")
            os.makedirs(nr_directory, exist_ok=True)
            print(nr_directory)
            command = ["update_blastdb.pl", "--decompress", self.blastnr_database]
            subprocess.run(command, cwd=nr_directory, check=True)

        if self.download_blastnt:
            print(
                "Downloading ",
                self.blastnt_database,
                "Blast database for Nucleotide Screening from NCBI \n"
            )
            nt_directory = os.path.join(self.database_directory,"nt_blast")
            os.makedirs(nt_directory, exist_ok=True)
            print(nr_directory)
            command = ["update_blastdb.pl", "--decompress", self.blastnt_database]
            subprocess.run(command, cwd=nt_directory, check=True)

        if self.download_taxonomy:
            print("Downloading Taxonomy databases from NCBI \n")
            tax_directory = os.path.join(self.database_directory,"taxonomy")
            os.makedirs(tax_directory, exist_ok=True)
            print(tax_directory)
            command = ["update_blastdb.pl", "--decompress", "taxdb"]
            subprocess.run(command, cwd=tax_directory, check=True)
            #command = ["wget","-c","-P",tax_directory, self.taxonomy_download_url]
            #result = subprocess.run(command, check = True)
            #if result.returncode != 0:
                #command_str = " ".join(command)
                #print(
                #    "\t ERROR: Command",
                #    command_str,
                #    "failed with error",
                #    result.stderr,
                #)

            # Parse the URL to extract the path
            #parsed_url = urllib.parse.urlparse(self.taxonomy_download_url)
            #taxonomy_filename_zipped = os.path.join(tax_directory,os.path.basename(parsed_url.path))

            #subprocess.run(["tar","-zxvf", taxonomy_filename_zipped], check=False)
            #os.remove(taxonomy_filename_zipped)

        print("\n\nThe common mechanism setup has completed!"
                "\nYou can find all downloaded databases in",
                self.database_directory,
                "\n\nHave a bio-safe and secure day!")

    def check_url_exists(self, url : str) -> bool:
        """ Helper function to quickly check if a URL is valid."""
        try:
            with urllib.request.urlopen(url) as response:
                # If the response status code is 200, the URL exists
                if response.status == 200:
                    return True
        except urllib.error.HTTPError as e:
            # Handle HTTP errors (like 404, 403, etc.)
            print(f"HTTP Error: {e.code}")
        except urllib.error.URLError as e:
            # Handle URL errors (like unreachable server, etc.)
            print(f"URL Error: {e.reason}")
        except ValueError as e:
            print("URL Value Error. It is likely the URL input is not a recognised URL format.")
        return False
    
    def print_step(self, i : int = 0, ii : int = -1):
        """ helper for quick step delinearation. """
        if ii > 0:
            print( "\n*----------------* Step ", i, ".", ii, "*----------------*")
            return
        print( "\n*----------------* Step ", i, "   *----------------*")
    
    def print_database_options(self):
        print("Fetching list of possible databases...")
        subprocess.run(["update_blastdb.pl", "--showall","pretty"], check=True)

    def user_input(self, prompt : str = ">>> "):
        """ Get input from the user, and do some basic string sanitation. """
        return input(prompt).strip().lower()

    def stop(self):
        """ Gracefully exit with a message to the user."""
        print("Exiting setup for The Common Mechanism.")
        sys.exit()

def add_args(parser_obj: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """
    parser_obj.add_argument(
        "-a",
        "--auto",
        dest="automated",
        default=False,
        action="store_true",
        help="Don't ask for user input, and use default options for everything.",
    )
    return parser_obj

def run(args):
    """ Run CLI with an parsed argument parser input."""
    CliSetup(args.automated)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)