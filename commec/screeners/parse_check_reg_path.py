import re
import os
import subprocess
import argparse

def parse_blast_output(input_text):
    # Initialize variables to store extracted information
    query = ""
    start = ""
    end = ""
    identity = ""
    sequences = ""
    regulated_species = ""
    regulated_taxids = ""
    non_regulated_species = ""
    non_regulated_taxids = ""

    # Use regex to extract information
    query_match = re.search(r'Input query file: (.+)', input_text)
    if query_match:
        query = query_match.group(1)

    coords_match = re.search(r'at bases (\d+) - (\d+)', input_text)
    if coords_match:
        start, end = coords_match.groups()

    seq_match = re.search(r'Best match to sequence\(s\) (.+?) at bases', input_text)
    if seq_match:
        sequences = seq_match.group(1)

    identity_match = re.search(r'Percent identity to query: (.+)', input_text)
    if identity_match:
        identity = identity_match.group(1)

    reg_species_match = re.search(r'Regulated species: (.+)', input_text)
    if reg_species_match:
        regulated_species = reg_species_match.group(1)

    reg_taxids_match = re.search(r'Regulated taxid\(s\): (.+)', input_text)
    if reg_taxids_match:
        regulated_taxids = reg_taxids_match.group(1)

    non_reg_species_match = re.search(r'Non-regulated species: (.+)', input_text)
    if non_reg_species_match:
        non_regulated_species = non_reg_species_match.group(1)

    non_reg_taxids_match = re.search(r'Non-regulated taxid\(s\): (.+)', input_text)
    if non_reg_taxids_match:
        non_regulated_taxids = non_reg_taxids_match.group(1)

    # Construct TSV line
    tsv_line = f"{query}\t{start}\t{end}\t{identity}\t{sequences}\t{regulated_species}\t{regulated_taxids}\t{non_regulated_species}\t{non_regulated_taxids}"

    return tsv_line

def process_file(input_file, script_path, database_path):
    cmd = [
        script_path,
        "-i", input_file,
        "-d", database_path,
        "-t", "8"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return parse_blast_output(result.stdout)

def main():
    parser = argparse.ArgumentParser(description="Process multiple nr.blastx files and generate a TSV output.")
    parser.add_argument("input_dir", help="Directory containing nr.blastx files")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("--script", required=True, help="Path to check_reg_path.py script")
    parser.add_argument("--db", required=True, help="Path to cm_databases directory")
    args = parser.parse_args()

    # Header for the TSV file
    header = "query\tstart\tend\tidentity\tsequence(s)\tregulated_species\tregulated_taxids\tnon-regulated_species\tnon-regulated_taxids"

    with open(args.output_file, 'w') as outfile:
        outfile.write(header + '\n')

        for filename in os.listdir(args.input_dir):
            if filename.endswith(".nr.blastx"):
                input_file = os.path.join(args.input_dir, filename)
                print("Processing input file {input_file}...")
                tsv_line = process_file(input_file, args.script, args.db)
                outfile.write(tsv_line + '\n')

    print(f"Processing complete. Output written to '{args.output_file}'")

if __name__ == "__main__":
    main()