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
        next(file)  # skip first line
        for line in file:
            fields = line.rstrip().split("\t")
            json_str += '    \\"cpu_time_seconds\\": ' + fields[0] + ',\n'
            json_str += '    \\"wall_clock_time\\": \\"' + fields[1] + '\\",\n'
            json_str += '    \\"mem_rss_mb\\": ' + fields[2] + ',\n'
            json_str += '    \\"mem_vms_mb\\": ' + fields[3] + ',\n'
            json_str += '    \\"mem_uss_mb\\": ' + fields[4] + ',\n'
            json_str += '    \\"mem_pss_mb\\": ' + fields[5] + ',\n'
            json_str += '    \\"io_in_bytes\\": ' + fields[6] + ',\n'
            json_str += '    \\"io_out_bytes\\": ' + fields[7] + ',\n'
            json_str += '    \\"mean_cpu_load\\": ' + fields[8] + ',\n'
            json_str += '    \\"cpu_load\\": ' + fields[9]
    return json_str
