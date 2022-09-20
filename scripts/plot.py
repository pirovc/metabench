#!/usr/bin/env python3
import argparse
import glob
import sys
import os
import json

import pandas as pd

# Bokeh
from bokeh.io import save
from bokeh.plotting import output_file
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.layouts import column


def main():

    parser = argparse.ArgumentParser(
        prog="plot metabench with bokeh", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input",     metavar="",
                        type=str, nargs="*", required=True, help="json file(s) and/or folder(s) with json reports (recursive search)")
    parser.add_argument("-o", "--output",    metavar="",
                        type=str, help="output html file")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s 1.0.0")

    if len(sys.argv) == 1:
        parser.print_help()
        return False

    args = parser.parse_args()

    json_files = find_json_files(args.input)
    parsed_json = parse_json(json_files)


    report_types = ["build", "binning", "profiling"]
    report_categories = ["benchmark", "size", "stats", "evals"]
    report_elements = ["config", "metrics"]
    tables = {}  # tables[type][category][element]
    for t in report_types:
        tables[t] = {}
        for c in report_categories:
            tables[t][c] = {}
            for e in report_elements:
                tables[t][c][e] = pd.DataFrame()


    for pjson in parsed_json.values():
        rep_type = pjson["report"]
        rep_cate = pjson["category"]
        tables[rep_type][rep_cate]["config"] =  pd.concat(
                    [tables[rep_type][rep_cate]["config"], load_config(pjson)], ignore_index=True)
        tables[rep_type][rep_cate]["metrics"] =  pd.concat(
                    [tables[rep_type][rep_cate]["metrics"], load_metrics(pjson)], ignore_index=True)
    

    for t in report_types:
        for c in report_categories:
            for e in report_elements:
                print(t,c,e)
                print(tables[t][c][e])

    

    # cds_config = ColumnDataSource(binning_config_evals_table.reindex())
    # columns = [
    #     TableColumn(field="tool", title="Tool"),
    #     TableColumn(field="version", title="Version"),
    #     TableColumn(field="sample", title="Sample"),
    #     TableColumn(field="database", title="DB"),
    #     TableColumn(field="database_arguments", title="DB Args"),
    #     TableColumn(field="arguments", title="Args"),
    # ]
    # data_table = DataTable(source=cds_config, columns=columns)

    # layout = column(data_table)
    # output_file(args.output, title="Metabench", mode="inline")
    # save(layout)

    return True


def find_json_files(input_list):
    files = []
    for i in input_list:
        if i.endswith(".json"):
            files.append(i)
        elif os.path.exists(i):
            for file in glob.glob(i + '/**/*.json', recursive=True):
                files.append(file)
    return files


def parse_json(json_files):
    parsed = {}
    for file in json_files:
        with open(file, "r") as jfile:
            parsed[file] = json.load(jfile)
    return parsed


def load_config(pjson):
    return pd.DataFrame([pjson["config"].values()], columns=pjson["config"].keys())

def load_metrics(pjson):
    cols = []
    vals = []
    for m1, v1 in pjson["metrics"].items():
        if isinstance(v1, dict):
            for m2, v2 in v1.items():
                if isinstance(v2, dict):
                    for m3, v3 in v2.items():
                        cols.append((m1, m2, m3))
                        vals.append(v3)
                else:
                    cols.append((m1, m2))
                    vals.append(v2)
        else:
            cols.append(m1)
            vals.append(v1)

    return pd.DataFrame([vals], columns=cols)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
