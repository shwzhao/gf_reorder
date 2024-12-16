#!/usr/bin/env python3

import argparse
from commands import match, reorder

def main():
    parser = argparse.ArgumentParser(prog='gf_reorder', description='GFF and FASTA files\' rename and reorder.')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    match.setup_parser(subparsers)
    reorder.setup_parser(subparsers)

    args = parser.parse_args()

    if args.subcommand:
        if args.subcommand == 'match':
            match.run(args)
        elif args.subcommand == 'reorder':
            reorder.run(args)
        else:
            parser.print_help()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
