workdir: config["workdir"]
include: "util.py"
include: "../tools/bracken.smk"
include: "../tools/ganon.smk"
include: "../tools/kmcp.smk"
include: "../tools/kraken2.smk"
include: "../tools/metacache.smk"

def input_build():
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                if vers in config["run"][tool]:
                    for dtbs in config["dbs"]:
                        if dtbs in config["run"][tool][vers]:
                            # define path to output
                            path = tool + "/" + vers + "/" + dtbs + "/"
                            # build product of all arguments (single or range)
                            for args in args_product(config["run"][tool][vers][dtbs]["args"] if "args" in config["run"][tool][vers][dtbs] else None):
                                out.append(path + join_args(sorted(args)))
    #import pprint
    #pprint.pprint(out)
    return out

rule all:
    input:
        build = expand("{i}.{ext}", i=input_build(), ext=["build.bench.json"])

rule bench:
    input:
        bench = "{tool}/{vers}/{dtbs}/{args}.build.bench.tsv",
        fsize = "{tool}/{vers}/{dtbs}/{args}.build.size.tsv"
    output:
        json = "{tool}/{vers}/{dtbs}/{args}.build.bench.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.args),
                                    "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers][wildcards.dtbs]["fixed_args"])}
    run:
        out_json = json_benchmark(input.bench, report="build", config=params.config)
        out_json["bench"]["total_size_bytes"] = 0
        with open(input.fsize, "r") as file:
            for line in file:
                s, f = line.rstrip().split("\t")
                out_json["bench"]["total_size_bytes"] += int(s)
        json_write(out_json, output.json)
