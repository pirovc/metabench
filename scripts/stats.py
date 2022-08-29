#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import json
from multitax import CustomTx, DummyTx, NcbiTx, GtdbTx


def main():

    default_ranks = ["superkingdom",
                     "phylum",
                     "class",
                     "order",
                     "family",
                     "genus",
                     "species"]

    parser = argparse.ArgumentParser(prog="stats", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input-file",     metavar="", type=file_exists, required=True, help="")
    parser.add_argument("-o", "--output-file",    metavar="", type=file_exists, help="json output file or STDOUT")
    parser.add_argument("-t", "--taxonomy",       metavar="", type=str, help="custom, ncbi or gtdb")
    parser.add_argument("-x", "--taxonomy-files", metavar="", type=file_exists, nargs="*", help="")
    parser.add_argument("-r", "--ranks",          metavar="", type=str, default=default_ranks, nargs="*", help="all for all ranks. empty for default ranks (superkingdom phylum class order family genus species)")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0.0")

    if len(sys.argv) == 1:
        parser.print_help()
        return False

    args = parser.parse_args()

    # Parse taxonomy if provided
    tax_args = {"undefined_node": "",
                "undefined_rank": "na",
                "undefined_name": "na",
                "root_name": "root",
                "root_rank": "root"}

    if args.taxonomy == "ncbi":
        tax = NcbiTx(files=args.taxonomy_files, **tax_args)
    elif args.taxonomy == "gtdb":
        tax = GtdbTx(files=args.taxonomy_files, **tax_args)
    elif args.taxonomy == "custom":
        tax = CustomTx(files=args.taxonomy_files,
                       cols=["node", "parent", "rank", "name"],
                       **tax_args)
    else:
        tax = DummyTx()

    # Generate results dict
    res = {"reads_classified": 0,
           "reads_invalid_tax": 0,
           "ranks": {}}

    # bioboxes:
    # headers start with @
    # fields: readid <tab> binid <tab> taxid
    for line in gzip.open(args.input_file, "rt") if args.input_file.endswith(".gz") else open(args.input_file,"r"):
        if line[0] == "@":
            continue
        fields = line.rstrip().split("\t")
        taxid = fields[2]

        if tax.latest(taxid) == tax.undefined_node:
            res["reads_invalid_tax"] += 1
            sys.stderr.write(taxid + " not found in taxonomy\n")
            continue
        elif taxid == tax.root_node:
            res["reads_invalid_tax"] += 1
            sys.stderr.write(taxid + " is root\n")
            continue
        else:
            rank = tax.rank(taxid)
            # if not in chosen ranks, get closest rank
            if args.ranks != ["all"] and rank not in args.ranks:
                # check linear reverse until find valid node
                for t in tax.lineage(taxid, ranks=args.ranks)[::-1]:
                    if t != tax.undefined_node:
                        rank = tax.rank(t)
                        break
                if rank not in args.ranks:
                    res["reads_invalid_tax"] += 1
                    sys.stderr.write(taxid + " without valid rank\n")
                    continue

            # account for rank
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


def file_exists(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError("{0} does not exist".format(file))
    return file


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
