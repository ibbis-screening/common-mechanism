#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science 
"""
Module for CLI setup of Commec, such that required 
databases are downloaded in a desired database directory.
"""
import sys
#import os
import argparse
import subprocess
import urllib.request

DESCRIPTION = "Helper script for downloading the databases required for running the Common Mechanism Screen"

class CliSetup:
    """
    Interacts with the user via the CLI, and then downloads appropriate files
    into appropraite folders as requested by the user, to allow for the Commec
    Screen workflow.
    """

    def __init__(self, automate : bool = False):
        self.ask : bool = not automate
        self.database_directory : str = "./commec-dbs/"
        self.download_biorisk : bool = True
        self.biorisk_download_url : str = "https://f005.backblazeb2.com/file/common-mechanism-dbs/common-mechanism-dbs.zip"
        self.download_blastnr : bool = False
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
              "\ninternet connection, and wget.",
              "\n -> You can exit this setup at any time with \"exit\"",
              "\n -> You can return to a previous step with \"back\"\n")
        
        self.setup_overall_directory()
        self.do_setup()
        
    def setup_overall_directory(self):
        """
        Get user inputs to ensure the
        """
        user_input : str = ""
        print("\nPlease provide the absolute or relative filepath",
              "to where you would like the Commec databases to be located...",
              "\nPress <Enter> to use default: ",
              self.database_directory)
        while(True):
            user_input : str = input("").strip()
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

            print("Using database directory: %s", self.database_directory)
            self.get_biorisk_url()
            return

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
        
        return False

    def get_biorisk_url(self):
        """
        Get the URL where the BIORISK database is located.
        """
        user_input : str = ""
        print("Please provide the URL to download the BIORISK database.",
              "\nPress <Enter> to use default: ",
              self.biorisk_download_url)
        while(True):
            user_input : str = input("").strip()
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
            
            print("Using Biorisk URL: %s", self.database_directory)
            self.confirm()
            return
        
    def confirm(self):
        """ Simply allows the user one last chance to confirm their settings."""
        print("The following settings will be used to setup Commec:",
              "\n -> Database Directory: ", self.database_directory,
              "\n -> BIORISK db URL: ", self.biorisk_download_url)
        print("\n\nPress <Enter> to confirm these settings, ",
              "\ntype \"back\" to alter previous settings, ",
              "\ntype \"exit\" to abort setup.",
              "\ntype \"restart\" to go back to the beginning.")
        while(True):
            user_input : str = input("").strip()
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

        if self.download_biorisk:
            command = ["wget", self.biorisk_download_url]
            print(
                "Downloading Biorisk database from \n%s", 
                self.biorisk_download_url
            )
            result = subprocess.run(command, check = True)
            if result.returncode != 0:
                command_str = " ".join(command)
                print(
                    "\t ERROR: Command %s failed with error %s",
                    command_str,
                    result.stderr,
                )

        if self.download_blastnr:
            print(
                "Downloading Biorisk database from \n%s", 
                self.biorisk_download_url
            )

        print("The common mechanism setup has completed!"
              " You can find all downloaded databases in %s",
                self.database_directory,
                "\nHave a bio-safe and secure day!")
        
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    CliSetup(args.automated)
