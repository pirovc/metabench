from itertools import product
import numpy as np
import datetime
import json

version = "1.0.0"
default_value = "default"
build_size_cmd = "du --bytes --dereference --max-depth 0"

def join_args(args):
    # Skip '' empty args created when using boolean parameters
    # If all are off, return default value
    j = "_".join([a for a in args if a != ""])
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
        return str(arg) if val == True else ""
    else:
        # If normal argument was deactivated --value ''
        if str(val) == "":
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
                expanded_args.append(
                    tuple([argval2str(p, element) for element in v]))
            elif isinstance(v, str):
                # Range values "--range-param": "0,1,0.1"
                if len(v.split(",")) == 3:
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


def json_default(report: str = "", config: dict = {}):
    j = {}
    j["report"] = report
    j["created"] = datetime.datetime.now().isoformat()
    j["config"] = config
    j["bench"] = {}
    return j


def json_benchmark(bench_file, report: str = "", config: dict = {}, sum_bench: dict = {}):
    out_json = json_default(report=report, config=config)
    out_json["config"] = config

    bench_fields = {0: "cpu_time_seconds",
                    # 1: "wall_clock_time",
                    2: "mem_rss_mb",
                    # 3: "mem_vms_mb",
                    # 4: "mem_uss_mb",
                    # 5: "mem_pss_mb",
                    6: "io_in_mb",
                    7: "io_out_mb",
                    8: "mean_cpu_load",
                    9: "cpu_load"}

    # Parse benchmark file and select fastest if more then one run
    with open(bench_file, "r") as file:
        next(file)  # skip first line with header
        bench_values = []
        for line in file:
            fields = line.rstrip().split("\t")
            bench_values.append(fields)

    # Select lowest execution time cpu_time_seconds
    selected_fields = bench_select_lowest(bench_values, 0)

    # Add repeats
    out_json["bench"]["repeats"] = len(bench_values)
    for field, metric in bench_fields.items():
        out_json["bench"][metric] = float(selected_fields[field])

    # If there's another benchmark to sum/max
    if sum_bench:
        for field, metric in bench_fields.items():
            if field == 0:  # sum cpu_time_seconds
                out_json["bench"][metric] += float(sum_bench["bench"][metric])
            else:  # max (any other field)
                out_json["bench"][metric] = max(
                    [out_json["bench"][metric], sum_bench["bench"][metric]])

    return out_json


def bench_select_lowest(bench_values, field):
    # No repetition
    if len(bench_values) == 1:
        return bench_values[field]
    else:
        # return min time
        val = [float(fields[field]) for fields in bench_values]
        min_idx = val.index(min(val))
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


