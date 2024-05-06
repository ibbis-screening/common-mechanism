"""
Command-line entrypoint for the package.
"""
import argparse
from common_mechanism.split import run as split_run, DESCRIPTION as split_DESCRIPTION, add_args as split_add_args

def main():
    parser = argparse.ArgumentParser(
        prog='commec',
        description="Command-line entrypoint for the Common Mechanism"
    )
    subparsers = parser.add_subparsers(dest='command')

    # Sub-command for "split"
    split_parser = subparsers.add_parser('split', help=split_DESCRIPTION)
    split_add_args(split_parser)

    args = parser.parse_args()
    if args.command == 'split':
        split_run(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
