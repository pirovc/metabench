#!/usr/bin/env python3
import argparse
import sys
import gzip
import json
from util import default_ranks, parse_tax, file_exists


def main():

    parser = argparse.ArgumentParser(
        prog="stats", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input-file",     metavar="",
                        type=file_exists, required=True, help="")
    parser.add_argument("-o", "--output-file",    metavar="",
                        type=file_exists, help="json output file or STDOUT")
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
    fixed_ranks = ["root"] + args.ranks if args.ranks else default_ranks

    # Generate results dict
    res = {"reads_classified": 0,
           "reads_invalid_tax": 0,
           "ranks": {}}

    # bioboxes:
    # headers start with @
    # fields: readid <tab> binid <tab> taxid
    for line in gzip.open(args.input_file, "rt") if args.input_file.endswith(".gz") else open(args.input_file, "r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = tax.latest(fields[2])

        if taxid == tax.undefined_node:
            res["reads_invalid_tax"] += 1
            sys.stderr.write(fields[2] + " not found in taxonomy\n")
            continue
        else:
            closest_taxid = tax.closest_parent(taxid, fixed_ranks)
            if closest_taxid == tax.undefined_node:
                res["reads_invalid_tax"] += 1
                sys.stderr.write(taxid + " without valid rank\n")
                continue

            # account for rank
            rank = tax.rank(closest_taxid)
            if rank not in res["ranks"]:
                res["ranks"][rank] = 0
            res["ranks"][rank] += 1
            res["reads_classified"] += 1

    if args.output_file:
        with open(args.output_file, "w") as outfile:
            json.dump(res, outfile, indent=4)
    else:
        print(json.dumps(res, indent=4))

    return True


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
