import os
import argparse
from multitax import CustomTx, DummyTx, NcbiTx, GtdbTx

default_ranks = ["root",
                 "superkingdom",
                 "phylum",
                 "class",
                 "order",
                 "family",
                 "genus",
                 "species"]


def parse_tax(tax, files):
    tax_args = {"undefined_node": "",
                "undefined_rank": "na",
                "undefined_name": "na",
                "root_name": "root",
                "root_rank": "root"}

    if tax == "ncbi":
        tax = NcbiTx(files=files, **tax_args)
    elif tax == "gtdb":
        tax = GtdbTx(files=files, **tax_args)
    elif tax == "custom":
        tax = CustomTx(files=files,
                       cols=["node", "parent", "rank", "name"],
                       **tax_args)
    else:
        tax = DummyTx()

    return tax


def file_exists(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError("{0} does not exist".format(file))
    return file
