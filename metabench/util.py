def join_params(params):
    return "_".join(params)


def params2str(arg, val):
    return str(arg) + "=" + str(val)


def str2params(s):
    return " ".join(s.replace("=", " ").split("_"))


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

def params_product(params):
    expanded_params = []
    for p, v in params.items():
        if isinstance(v, list):
            rparams = []
            # include stop
            for range_val in list(np.arange(v[0], v[1], v[2])) + [v[1]]:
                # round to 2 decimal
                rparams.append(params2str(p, round(range_val, 2)))
            expanded_params.append(tuple(rparams))
        else:
            expanded_params.append(tuple([params2str(p, v)]))

    prod_params = []
    if expanded_params:
        prod_params = list(product(*expanded_params))
    else:
        prod_params = ["default"]

    return prod_params
