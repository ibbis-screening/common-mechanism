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
import ftplib
from urllib import request, error, parse
import zipfile
import tarfile

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
        self.database_directory : str = "commec-dbs/"

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
        self.default_taxonomy_download_url : str = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
        self.taxonomy_download_url : str = self.default_taxonomy_download_url

        if automate:
            self.download_biorisk = True
            self.download_blastnr = True
            self.download_blastnt = True
            self.download_taxonomy = True
            self.do_setup()
        else:
            self.start()

    def start(self):
        """
        Starts the user interrogation process.
        """
        print("""\n                       Welcome to\n
 ██████╗ ██████╗ ███╗   ███╗███╗   ███╗███████╗ ██████╗ \033[38;5;202m         ▄▄               \033[0m
██╔════╝██╔═══██╗████╗ ████║████╗ ████║██╔════╝██╔════╝ \033[38;5;202m       ▄███▌              \033[0m
██║     ██║   ██║██╔████╔██║██╔████╔██║█████╗  ██║      \033[38;5;202m      ▐█████              \033[0m
██║     ██║   ██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██║      \033[38;5;202m     ▐██████▌             \033[0m
╚██████╗╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗╚██████╗ \033[38;5;202m     ███████▌             \033[0m
 ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝ ╚═════╝ \033[38;5;202m    ▐███████▌             \033[0m
\033[48;5;17m█ █████▄ █████▄ █ ▄█▀█▄                                 \033[38;5;202m     ███████   ▄█▄      \033[0m
\033[48;5;17m█ █    █ █    █ █ █   ▀          DATABASE               \033[38;5;202m      █████▌  ▄███▄▄     \033[0m  
\033[48;5;17m█ █████▄ █████▄ █ ▀███▄            SETUP                \033[38;5;202m      ▐█████▄██▀    ▀▄    \033[0m   
\033[48;5;17m█ █    █ █    █ █ ▄   █              UTILITY            \033[38;5;202m      ▐████████       ▌  \033[0m
\033[48;5;17m█ █████▀ █████▀ █ ▀█▄█▀                                 \033[38;5;202m      ████████▀         \033[0m
                                                        \033[38;5;202m   ▄▄██████▀▀             \033[0m
                                                        \033[38;5;202m ▀▀                       \033[0m""")
        print("                 The Common Mechanism!",
              "\n\nInternational Biosecurity and Biosafety Initiative for Science",
              "\nCopyright © 2021-2024 ",
              "\n\nThis script will help download the mandatory databases ",
              "\nrequired for using Commec Screen, and requires a stable",
              "\ninternet connection, wget, and update_blastdb.pl.")
        self.print_help_info()
        print()
        self.check_requirements()
        self.setup_overall_directory()
        self.do_setup()

    def print_help_info(self, additional_help = [""]):
        """ 
        Prints helpful instructions according to the current step, 
        as well as general instructions. 
        """
        add_help_str = "" if len(additional_help) == 0 else "".join(additional_help)
        print("\033[38;5;202m"
            " Instructions: "
            "\n -> You can exit this setup at any time with \"exit\""
            "\n -> You can return to a previous step with \"back\""
            "\n -> You can get additional help with \"help\""
              + add_help_str +
            "\033[0m")

    def check_requirements(self):
        """ Checks for wget, and update_blastdb.pl, which should both be present
        given the conda environment to install and use Commec. """
        print("\033[38;5;242mChecking for wget, and update_blastdb")
        dont_proceed = False
        result = subprocess.run(["which", "wget"], check=False)
        if result.returncode != 0:
            print("\033[38;5;202mDependency wget could not be found. "
                  "Check wget is installed in this environment.")
            dont_proceed = True
        result = subprocess.run(["which", "update_blastdb.pl"], check=False)
        if result.returncode != 0:
            print("\033[38;5;202mDependency update_blastdb.pl could not be found. "
                  "Blast is likely not installed in this python environment.")
            dont_proceed = True
        if dont_proceed:
            print("\033[38;5;202mError: Required dependency missing!")
            self.stop()
        print("\033[0m")

    def check_directory_is_writable(self, input_directory : str) -> bool:
        """ Checks a directory is viable by creating it and destroying it."""
        try:
            os.makedirs(input_directory, exist_ok=True)
            try: # Removedirs can fail if it is a single leaf, containing other files.
                os.removedirs(input_directory)
            except OSError:
                pass
        except subprocess.CalledProcessError:
            return False
        return True

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
        while True:
            user_input : str = self.user_input()
            if user_input == "back":
                self.stop() # This is the first step, no point going back.
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> Provide a path to download into e.g. \"/location/to/download/commec-dbs\"",
                                      "\n (The full path will be created if it does not exist.)"])
                continue

            if len(user_input) > 0:
                directory_is_valid = self.check_directory_is_writable(user_input)
                if not directory_is_valid:
                    print(user_input, " is not a valid directory structure!")
                    continue
                self.database_directory = user_input

            print("Using database directory: ", self.database_directory)
            self.decide_commec_dbs()
            return

    def decide_commec_dbs(self):
        """ Decide whether the Commec Benign/risks database needs to be downloaded. """
        self.print_step(2)
        print("Do you want to download the "
              "mandatory Commec databases? (~1.2 GB)",
              "\n\"y\" or \"n\", for yes or no.")
        while True:
            user_input : str = self.user_input()
            if user_input == "help":
                self.print_help_info(["\n -> type yes,y or no,n to indicate decision.",
                                      "\n  (The Commec databases consist of a currated biorisk",
                                      "\n  and benign database, which are required for commec to run",
                                      "\n  and are the only databases used in \"--fast-mode\")"])
                continue
            if user_input == "back":
                self.setup_overall_directory()
                return
            if user_input in ["y","yes"]:
                self.download_biorisk = True
                self.get_biorisk_url()
                return
            if user_input in ["n","no"]:
                self.download_biorisk = False
                self.decide_blastnr()
                return
            print("Unrecognised input (", user_input, ")")

    def get_biorisk_url(self):
        """
        Get the URL where the Commec Biorisk and Benign databases are located.
        """
        self.print_step(2,1)
        print("Please provide the URL to download the Commec database.",
              "\nPress <Enter> to use existing: ",
              self.biorisk_download_url)
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> Provide a URL to download, the commec database URL can be found at ",
                                      "\n  https://github.com/ibbis-screening/common-mechanism/wiki/install",
                                      "\n  (The URL will be checked that it is valid",
                                      "\n  which will require an internet connection.)"])
                continue
            if user_input == "back":
                self.decide_commec_dbs()
                return
            if len(user_input) > 0:
                self.biorisk_download_url = user_input

            print("Checking URL is valid ... ")
            if not self.check_url_exists(self.biorisk_download_url):
                print(self.biorisk_download_url, " is not a valid URL! "
                      "(or you are not connected to the internet)")
                continue

            print("Using Commec Biorisk and Benign URL:", self.biorisk_download_url)
            self.decide_blastnr()
            return

    def decide_blastnr(self):
        """ Decide whether a Protein database needs to be downloaded. """
        self.print_step(3)
        print("Do you want to download the"
              " protein NR database for protein screening? (~530 GB)",
              "\n\"y\" or \"n\", for yes or no.")
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> type yes,y or no,n to indicate decision.",
                                      "\n   (The Blast non-redundant protein database is used to screen",
                                      "\n   translated queries when not in \"--fast-mode\", and is ",
                                      "\n   around 530 GB in size.)"])
                continue
            if user_input == "back":
                if self.download_biorisk:
                    self.get_biorisk_url()
                    return
                self.decide_commec_dbs()
                return
            if user_input == "y" or user_input == "yes":
                self.download_blastnr = True
                self.decide_blastnt()
                return
            if user_input == "n" or user_input == "no":
                self.download_blastnr = False
                self.decide_blastnt()
                return
            print("Unrecognised input (", user_input, ")")

    def decide_blastnt(self):
        """ Decide what Nucleotide database needs to be downloaded. """
        self.print_step(4)
        print("Do you want to download the Nucleotide NT databases for non-coding"
              " region nucleotide screening? (~580 GB)",
              "\n\"y\" or \"n\", for yes or no.")
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> type yes,y or no,n to indicate decision.",
                                      "\n   (The Blast Nucleotide database is used to screen",
                                      "\n   non-coding regions of queries when not in \"--fast-mode\","
                                      "\n   and is around 580 GB in size.)"])
                continue
            if user_input == "back":
                self.decide_blastnr()
                return
            if user_input == "y" or user_input == "yes":
                self.download_blastnt = True
                self.decide_taxonomy()
                return
            if user_input == "n" or user_input == "no":
                self.download_blastnt = False
                self.decide_taxonomy()
                return
            print("Unrecognised input (", user_input, ")")

    def decide_taxonomy(self):
        """ Decide whether taxonomy database need to be downloaded. """
        self.print_step(5)
        print("Do you want to download the Taxonomy databases? ( less than ~500 MB)",
              "\n\"y\" or \"n\", for yes or no.")
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> type yes,y or no,n to indicate decision.",
                                      "\n   (A Taxonomy database is used to check taxonomy IDs",
                                      "\n   for regulation when not in \"--fast-mode\","
                                      "\n   and is around 500 MB in size.)"])
                continue
            if user_input == "back":
                self.decide_blastnt()
                return
            if user_input == "y" or user_input == "yes":
                self.download_taxonomy = True
                self.get_taxonomy_url()
                return
            if user_input == "n" or user_input == "no":
                self.download_taxonomy = False
                self.confirm()
                return
            print("Unrecognised input (", user_input, ")")

    def get_taxonomy_url(self):
        """
        Get the URL where the Commec Biorisk and Benign databases are located.
        """
        self.print_step(5,1)
        user_input : str = ""
        print("Please provide the URL to download the Taxonomy database.",
              "\nPress <Enter> to use existing: ",
              self.taxonomy_download_url)
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> Provide a URL to download, the taxonomy URL can be found at ",
                                      "\n  https://github.com/ibbis-screening/common-mechanism/wiki/install",
                                      "\n  (The URL will be checked that it is valid",
                                      "\n  which will require an internet connection.)"])
                continue
            if user_input == "back":
                self.decide_taxonomy()
                return
            if len(user_input) > 0:
                self.taxonomy_download_url = user_input
                continue

            print("Checking URL is valid ... ")
            if not self.check_url_exists(self.taxonomy_download_url):
                print(self.taxonomy_download_url, " is not a valid URL! "
                      "(or you are not connected to the internet)")
                continue

            print("Using Commec Biorisk and Benign URL:", self.taxonomy_download_url)
            self.confirm()
            return

    def confirm(self):
        """ Simply allows the user one last chance to confirm their settings."""
        self.print_step(6)
        print("The following settings will be used to setup Commec:",
              "\n -> Database Directory: ", self.database_directory)
        if self.download_biorisk:
            print(" -> Commec Biorisk and Benign databases will be downloaded,"
                  "\n    from URL: ", self.biorisk_download_url)
        if self.download_blastnr:
            print(" -> Protein NR database will be downloaded.")
        if self.download_blastnt:
            print(" -> Nucleotide NT database will be downloaded.")
        if self.download_taxonomy:
            print(" -> Taxonomy database will be downloaded,"
                  "\n    from URL: ", self.taxonomy_download_url)

        if not (self.download_biorisk or
                self.download_blastnr or
                self.download_blastnt or
                self.download_taxonomy):
            print("You have opted to not download anything!"
                  "\nContinuing will simply exit this setup. "
                  "\nGo back to confirm a database to download.")

        print("\nPress <Enter> to confirm these settings, ",
              "\ntype \"back\" to alter previous setting, ",
              "\ntype \"exit\" to abort setup.",
              "\ntype \"start\" to alter from the beginning.")
        while True:
            user_input : str = self.user_input()
            if user_input in ["help", "h"]:
                self.print_help_info(["\n -> press <Enter> to confirm choices, and start the download process.",
                                      "\n -> type \"start\" to alter choices from the beginning."])
                continue
            if len(user_input) == 0:
                return
            if user_input == "back":
                if self.download_taxonomy:
                    self.get_taxonomy_url()
                    return
                self.decide_taxonomy()
                return
            if user_input == "start":
                self.setup_overall_directory()
                return
            print("Unrecognised input (", user_input, ")")

    def do_setup(self):
        """ 
        Once CliSetup state has been finalized,
        call to perform the required actions. 
        """
        self.print_step(7)
        if not (self.download_biorisk or
                self.download_blastnr or
                self.download_blastnt or
                self.download_taxonomy):
            print("No downloads were requested!")
            self.stop()

        os.makedirs(self.database_directory, exist_ok=True)

        if self.download_biorisk:
            command = [
                "wget","-c","-P",
                self.database_directory,
                self.biorisk_download_url
            ]
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
            parsed_url = parse.urlparse(self.biorisk_download_url)
            filename_zipped = os.path.join(
                self.database_directory,
                os.path.basename(parsed_url.path)
            )

            print("Extracting Biorisk databases...")
            # Open the zip file and extract its contents, remove zip file.
            with zipfile.ZipFile(filename_zipped, 'r') as zip_ref:
                zip_ref.extractall(self.database_directory)
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

            command = [
                "wget", "-c", "-P", tax_directory, self.taxonomy_download_url
            ]
            print("Downloading taxonomy data from", self.taxonomy_download_url)
            result = subprocess.run(command, check=True)

            # Check if the download was successful
            if result.returncode != 0:
                command_str = " ".join(command)
                print(f"\t ERROR: Command {command_str} failed with error {result.stderr}")

            # Parse the URL to get the filename
            parsed_url = parse.urlparse(self.taxonomy_download_url)
            filename_zipped = os.path.join(tax_directory, 
                                           os.path.basename(parsed_url.path))

            # Extract the tar.gz file
            print("Extracting taxonomy databases...")
            with tarfile.open(filename_zipped, "r:gz") as tar:
                tar.extractall(path=tax_directory)

            os.remove(filename_zipped)

        print("\n\nThe common mechanism setup has completed!"
                "\nYou can find all downloaded databases in",
                self.database_directory,
                "\n\nHave a bio-safe and secure day!")

    def check_url_exists(self, url : str) -> bool:
        """ Helper function to quickly check if a URL is valid."""
        parsed_url = parse.urlparse(url)

        if parsed_url.scheme in ['http', 'https']:
            # Handle HTTP/HTTPS URLs using urllib
            try:
                with request.urlopen(url) as response:
                    # If the response status code is 200, the URL exists
                    if response.status == 200:
                        return True
            except error.HTTPError as e:
                # Handle HTTP errors (like 404, 403, etc.)
                print(f"HTTP Error: {e.code}")
            except error.URLError as e:
                # Handle URL errors (like unreachable server, etc.)
                print(f"URL Error: {e.reason}")
            except ValueError as e:
                print("URL Value Error. It is likely the URL input"
                    " is not a recognized URL format.")

        elif parsed_url.scheme == 'ftp':
            # Handle FTP URLs using ftplib
            try:
                with ftplib.FTP(parsed_url.hostname) as ftp:
                    ftp.login()
                    # Try to change to the directory and check the file
                    ftp.cwd(os.path.dirname(parsed_url.path))
                    ftp.size(os.path.basename(parsed_url.path))
                    return True
            except ftplib.error_perm as e:
                print(f"FTP Permission Error: {e}")
            except ftplib.error_temp as e:
                print(f"FTP Temporary Error: {e}")
            except ftplib.all_errors as e:
                print(f"FTP Error: {str(e)}")

        return False


    def print_step(self, i : int = 0, ii : int = -1):
        """ helper for quick step delinearation. """
        if ii > 0:
            print( "\n\033[38;5;17m*----------------*\033[0m Step ",
                   i, ".", ii, "\033[38;5;17m*----------------*\033[0m")
            return
        print( "\n\033[38;5;17m*----------------*\033[0m Step ",
               i, " \033[38;5;17m*-------------------*\033[0m")

    def print_database_options(self):
        """ 
        Fetches the databases from NCBI that can be downloaded,
        as well as their sizes, and descriptions. 
        """
        print("Fetching list of possible databases...")
        subprocess.run(["update_blastdb.pl", "--showall","pretty"], check=True)

    def user_input(self, prompt : str = ">>> "):
        """ Get input from the user, and do some basic string sanitation. """
        user_input : str = input(prompt).strip().lower()
        if user_input == "exit":
            self.stop()
        return user_input

    def stop(self):
        """ Gracefully exit with a message to the user."""
        print("\033[0mExiting setup for The Common Mechanism.")
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

def run(arguments):
    """ Run CLI with an parsed argument parser input."""
    CliSetup(arguments.automated)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)
