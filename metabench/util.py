from itertools import product
import numpy as np


def join_args(args):
    return "_".join(args)


def dict2args(args_dict):
    s = ""
    if args_dict is not None:
        for arg, val in args_dict.items():
            s += str(arg) + " " + str(val) + " "
    return s


def str2args(s):
    if s == "default":
        return ""
    else:
        return " ".join(s.replace("=", " ").split("_"))


def argval2str(arg, val):
    return str(arg) + "=" + str(val)


def args_product(config_args_dict):
    prod_args = [["default"]]
    if config_args_dict is not None:
        expanded_args = []
        for p, v in config_args_dict.items():
            if isinstance(v, list):
                rargs = []
                # include stop
                for range_val in list(np.arange(v[0], v[1], v[2])) + [v[1]]:
                    # round to 2 decimal
                    rargs.append(argval2str(p, round(range_val, 2)))
                expanded_args.append(tuple(rargs))
            else:
                expanded_args.append(tuple([argval2str(p, v)]))

        if expanded_args:
            prod_args = list(product(*expanded_args))

    return prod_args


def json_wildcards(wdict):
    json_str = ""
    for p, v in wdict.items():
        json_str += '\\"' + p + '\\": \\"' + v + '\\",\n'
    return json_str


def json_bench(bench_file):
    json_str = ""
    with open(bench_file, "r") as file:
        next(file)  # skip first line with header

        bench_values = []
        for line in file:
            fields = line.rstrip().split("\t")
            bench_values.append(fields)

    selected_fields = bench_select(bench_values)

    json_str += '    \\"cpu_time_seconds\\": ' + selected_fields[0] + ',\n'
    json_str += '    \\"wall_clock_time\\": \\"' + selected_fields[1] + '\\",\n'
    json_str += '    \\"mem_rss_mb\\": ' + selected_fields[2] + ',\n'
    json_str += '    \\"mem_vms_mb\\": ' + selected_fields[3] + ',\n'
    json_str += '    \\"mem_uss_mb\\": ' + selected_fields[4] + ',\n'
    json_str += '    \\"mem_pss_mb\\": ' + selected_fields[5] + ',\n'
    json_str += '    \\"io_in_bytes\\": ' + selected_fields[6] + ',\n'
    json_str += '    \\"io_out_bytes\\": ' + selected_fields[7] + ',\n'
    json_str += '    \\"mean_cpu_load\\": ' + selected_fields[8] + ',\n'
    json_str += '    \\"cpu_load\\": ' + selected_fields[9]

    return json_str

def bench_select(bench_values):
    # No repetition
    if len(bench_values)==1:
        return bench_values[0]
    else:
        # return min time
        cpu_time_seconds = [float(fields[0]) for fields in bench_values]
        min_idx = cpu_time_seconds.index(min(cpu_time_seconds))
        return bench_values[min_idx]

def header_bioboxes_classify(tool, wildcards):
    #https://github.com/bioboxes/rfc/blob/master/data-format/binning.mkd
    header = ""
    header += "@Version:0.10.0\n"
    header += "@SampleID: " + tool + " " + " ".join(wildcards) + "\n"
    header += "@@SEQUENCEID\tTAXID\tBINID"
    return header


def header_bioboxes_profile(tool, wildcards):
    #https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd
    header = ""
    header += "@Version:0.10.0\n"
    header += "@SampleID: " + tool + " " + " ".join(wildcards) + "\n"
    header += "@Ranks: " + "\n"
    header += "@Taxonomy: " + "\n"
    header += "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
    return header