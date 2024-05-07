#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Command-line entrypoint for the package. Calls `screen.py`, `flag.py` and `split.py` as subcommands.

Command-line usage:
    - commec screen -d /path/to/databases input.fasta
    - commec flag /path/to/directory/with/output.screen 
    - commec split input.fasta
"""
import argparse
from common_mechanism.flag import (
    DESCRIPTION as flag_DESCRIPTION,
    add_args as flag_add_args,
    run as flag_run
)
from common_mechanism.split import (
    DESCRIPTION as split_DESCRIPTION,
    add_args as split_add_args,
    run as split_run
)

def main():
    """
    Parse the command line arguments and call the relevant sub-command.
    """
    parser = argparse.ArgumentParser(
        prog='commec',
        description="Command-line entrypoint for the Common Mechanism"
    )
    subparsers = parser.add_subparsers(dest='command')

    # Sub-command for "split"
    split_parser = subparsers.add_parser('split', help=split_DESCRIPTION)
    split_add_args(split_parser)

    # Sub-command for "flag"
    flag_parser = subparsers.add_parser('flag', help=flag_DESCRIPTION)
    flag_add_args(flag_parser)

    args = parser.parse_args()
    if args.command == 'split':
        split_run(args)
    elif args.command == 'flag':
        flag_run(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
