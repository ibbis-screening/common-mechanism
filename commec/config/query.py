"""
Container class to hold information pertaining to query from an input fasta file,
 as well as derived information, such as translated sequences, whether or not 
 the query was derived from AA or NT.
"""
import os
import subprocess

class Query:
    """
    Query to screen, based on an input FASTA. Self-calculates AA version.
    Future direction to back translate to NT when given AA too.
    """
    def __init__(self, input_fasta_filepath : str):
        self.input_fasta_path = input_fasta_filepath
        self.nt_path : str = ""
        self.aa_path : str = ""
        self.querie_names = []

    def validate(self, output_prefix : str):
        """ 
        Translate or reverse translate query, to be ready in AA or NT format. 
        """
        self.nt_path = Query.get_cleaned_fasta(self.input_fasta_path, output_prefix)
        self.aa_path = f"{output_prefix}.transeq.faa"
        self.get_query_names()

        # Run command transeq, to translate our input sequences.
        command = ["transeq", self.nt_path, self.aa_path, "-frame", "6", "-clean"]
        try:
            result = subprocess.run(command, check=True)
            if result.returncode != 0:
                command_str = ' '.join(command)
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                )
        except RuntimeError as error:
            raise RuntimeError("Input FASTA {fasta_to_screen} could not be translated.") from error
        
    def get_query_names(self):
        """ Populates query names with a list of all the queries within a fasta, for quick reference."""
        self.querie_names = []
        with (open(self.nt_path, "r", encoding="utf-8") as fin):
            for line in fin:
                if line.startswith('>'):
                    self.querie_names.append(line[1:-1])
        return
        
    @staticmethod
    def get_cleaned_fasta(input_file, out_prefix):
        """
        Return a FASTA where whitespace (including non-breaking spaces) and 
        illegal characters are replaced with underscores.
        """
        cleaned_file = f"{out_prefix}.cleaned.fasta"
        with (
            open(input_file, "r", encoding="utf-8") as fin,
            open(cleaned_file, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                modified_line = "".join(
                    "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c for c in line
                )
                fout.write(f"{modified_line}{os.linesep}")
        return cleaned_file