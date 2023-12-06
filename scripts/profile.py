#!/usr/bin/env python3
import argparse
import sys
import gzip
import json
from util import default_ranks, parse_tax, file_exists
from collections import OrderedDict


def main():

    parser = argparse.ArgumentParser(
        prog="stats", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input-file",     metavar="",
                        type=file_exists, required=True, help="")
    parser.add_argument("-x", "--taxid-col",    metavar="", type=int)
    parser.add_argument("-p", "--perc-col",    metavar="", type=int)
    parser.add_argument("-t", "--taxonomy",       metavar="",
                        type=str, help="custom, ncbi or gtdb")
    parser.add_argument("-x", "--taxonomy-files", metavar="",
                        type=file_exists, nargs="*", help="")
    parser.add_argument("-r", "--ranks",          metavar="", type=str, nargs="*",
                        help="empty for default ranks (superkingdom phylum class order family genus species)")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s 1.0.0")

    if len(sys.argv) == 1:
        parser.print_help()
        return False

    args = parser.parse_args()

    tax = parse_tax(args.taxonomy, args.taxonomy_files)
    # without root
    fixed_ranks = args.ranks if args.ranks else default_ranks[1:]
    tax.build_lineages(ranks=fixed_ranks)

    with open(args.input_file, "r") as file:
        for line in file:
            fields = line.rstrip().split("\t")
            taxid = tax.latest(fields[args.taxid_col])
            rank = tax.rank(taxid)
            if taxid == tax.undefined_node or rank == tax.undefined_rank:
                sys.stderr.write(
                    "[" + line.rstrip() + "] skipping taxid/rank not defined in taxonomy\n")
            elif rank not in fixed_ranks:
                sys.stderr.write(
                    "[" + line.rstrip() + "] skipping - not in defined ranks\n")
            else:
                prof = ""
                prof += taxid + "\t"
                prof += rank + "\t"
                prof += "|".join(tax.lineage(taxid,
                                             ranks=fixed_ranks[0:fixed_ranks.index(rank)+1])) + "\t"
                prof += "|".join(tax.name_lineage(taxid,
                                                  ranks=fixed_ranks[0:fixed_ranks.index(rank)+1])) + "\t"
                prof += fields[args.perc_col] + "\n"
                sys.stdout.write(prof)

    return True


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
