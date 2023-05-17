workdir: config["workdir"]
include: "util.py"
include: "../tools/bracken.smk"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_all():
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                if vers in config["run"][tool]:
                    for samp in config["samples"]:
                        for dtbs in config["run"][tool][vers]["dbs"]:
                            # Get all folders in database prefix
                            db_dir = config["run"][tool][vers]["dbs"][dtbs]
                            for dtbs_args in [d for d in os.listdir(db_dir) if os.path.isdir(os.path.join(db_dir, d))]:
                                # define path to output
                                path = tool + "/" + vers + "/" + samp + "/" + dtbs + "/" + dtbs_args + "/"
                                # build product of all arguments (single or range)
                                for binning_args in args_product(config["run"][tool][vers]["binning_args"]):
                                    ba = join_args(sorted(binning_args))
                                    out.append(path + ba + ".binning")
                                    for profiling_args in args_product(config["run"][tool][vers]["profiling_args"]):
                                        pa = join_args(sorted(profiling_args))
                                        out.append(path + ba + "/" + pa + ".profiling")
    #import pprint
    #pprint.pprint(out)
    return out

rule all:
    input:
        expand("{i}.{ext}", i=input_all(), ext=["bioboxes", "bench.json"])

rule bench_binning:
    input:
        bench = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args),
                                    "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers]["fixed_args"])}
    run:
        json_write(json_benchmark(input.bench, report="binning", config = params.config), output.json)

rule bench_profiling:
    input:
        bench = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args),
                                    "profiling_arguments": str2args(wildcards.p_args),
                                    "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers]["fixed_args"])}
    run:
        json_write(json_benchmark(input.bench, report="profiling", config = params.config), output.json)