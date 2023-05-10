from itertools import product
import numpy as np
import datetime
import json

default_value = "default"

def join_args(args):
    # Skip '' empty args created when using boolean parameters
    # If all are off, return default value
    j = "_".join([a for a in args if a!=""])
    return j if j else default_value

def dict2args(args_dict):
    s = ""
    if args_dict is not None:
        for arg, val in args_dict.items():
            s += str(arg) + " " + str(val) + " "
    return s


def str2args(s):
    if s == default_value:
        return ""
    else:
        return " ".join(s.replace("=", " ").split("_"))


def argval2str(arg, val):
    if isinstance(val, bool):
        # Boolean argument are either single or none (--active or "")
        return str(arg) if val==True else ""
    else:
        # If normal argument was deactivated --value ''
        if str(val)=="":
            return ""
        else:
            return str(arg) + "=" + str(val)


def args_product(config_args_dict):
    prod_args = [[default_value]]
    if config_args_dict is not None:
        expanded_args = []
        for p, v in config_args_dict.items():
            # List of values
            # "--param": [0, 0.2, 0.5]
            # "--bool-param": [True, False]
            if isinstance(v, list):
                expanded_args.append(tuple([argval2str(p, element) for element in v]))
            elif isinstance(v, str):
                # Range values "--range-param": "0,1,0.1"
                if len(v.split(","))==3:
                    split_v = v.split(",")
                    rargs = []
                    # include stop
                    for range_val in list(np.arange(float(split_v[0]), float(split_v[1]), float(split_v[2]))) + [float(split_v[1])]:
                        # round to 2 decimal
                        rargs.append(argval2str(p, round(range_val, 2)))
                    expanded_args.append(tuple(rargs))
                else:
                    # Single string param --mode avg
                    expanded_args.append(tuple([argval2str(p, v)]))
            else:
                # Single value --cutoff 0.5
                expanded_args.append(tuple([argval2str(p, v)]))

        if expanded_args:
            prod_args = list(product(*expanded_args))

    return prod_args


def json_load(json_file):
    with open(json_file, "r") as file:
        return json.load(file)


def json_write(j, json_file):
    with open(json_file, "w") as file:
        json.dump(j, file, indent=4)


def json_default(report: str = "", category: str = "", config: dict = {}):
    j = {}
    j["report"] = report
    j["category"] = category
    j["created"] = datetime.datetime.now().isoformat()
    j["config"] = config
    j["metrics"] = {}
    return j


def json_benchmark(bench_file, report: str = "", category: str = "", config: dict = {}):
    out_json = json_default(report=report, category=category, config=config)
    out_json["config"] = config

    # Parse benchmark file and select fastest if more then one run
    with open(bench_file, "r") as file:
        next(file)  # skip first line with header
        bench_values = []
        for line in file:
            fields = line.rstrip().split("\t")
            bench_values.append(fields)

    selected_fields = bench_select(bench_values)
    out_json["metrics"]["repeats"] = len(bench_values)
    out_json["metrics"]["cpu_time_seconds"] = float(
        selected_fields[0])
    out_json["metrics"]["wall_clock_time"] = selected_fields[1]
    out_json["metrics"]["mem_rss_mb"] = float(selected_fields[2])
    out_json["metrics"]["mem_vms_mb"] = float(selected_fields[3])
    out_json["metrics"]["mem_uss_mb"] = float(selected_fields[4])
    out_json["metrics"]["mem_pss_mb"] = float(selected_fields[5])
    out_json["metrics"]["io_in_bytes"] = float(selected_fields[6])
    out_json["metrics"]["io_out_bytes"] = float(
        selected_fields[7])
    out_json["metrics"]["mean_cpu_load"] = float(
        selected_fields[8])
    out_json["metrics"]["cpu_load"] = float(selected_fields[9])
    return out_json


def bench_select(bench_values):
    # No repetition
    if len(bench_values) == 1:
        return bench_values[0]
    else:
        # return min time
        cpu_time_seconds = [float(fields[0]) for fields in bench_values]
        min_idx = cpu_time_seconds.index(min(cpu_time_seconds))
        return bench_values[min_idx]


def header_bioboxes_binning(tool, wildcards):
    # https://github.com/bioboxes/rfc/blob/master/data-format/binning.mkd
    header = ""
    header += "@Version:0.10.0\n"
    header += "@SampleID:" + tool + " " + " ".join(wildcards) + "\n"
    header += "@@SEQUENCEID\tTAXID\t__NCBI_ASSEMBLY_ACCESSION\t__NCBI_SEQUENCE_ACCESSION"
    return header


def header_bioboxes_profiling(tool, ranks, taxonomy_files, wildcards):
    # https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd
    header = ""
    header += "@Version:0.10.0\n"
    header += "@SampleID: " + tool + " " + " ".join(wildcards) + "\n"
    header += "@Ranks: " + "|".join(ranks) + "\n"
    header += "@Taxonomy: " + ",".join(taxonomy_files) + "\n"
    header += "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
    return header
