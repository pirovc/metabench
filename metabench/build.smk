workdir: config["workdir"]
include: "util.py"
include: "../tools/bracken.smk"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_build():
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                if vers in config["run"][tool]:
                    for dtbs in config["dbs"]:
                        # define path to output
                        path = tool + "/" + vers + "/" + dtbs + "/"
                        # build product of all arguments (single or range)
                        for args in args_product(config["run"][tool][vers][dtbs]["args"] if "args" in config["run"][tool][vers][dtbs] else None):
                            out.append(path + join_args(args))
    return out

rule all:
    input:
        build = expand("{i}.{ext}", i=input_build(), ext=["build.bench.json", "build.size.json"])

rule bench:
    input:
        bench = "{tool}/{vers}/{dtbs}/{args}.build.bench.tsv"
    output:
        json = "{tool}/{vers}/{dtbs}/{args}.build.bench.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "database": wildcards.dtbs,
                                    "arguments": str2args(wildcards.args),
                                    "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers][wildcards.dtbs]["fixed_args"])}
    run:
        json_write(json_benchmark(input.bench, mode="build", category="benchmark", config=params.config), output.json)

rule size:
    input:
        fsize = "{tool}/{vers}/{dtbs}/{args}.build.size.tsv"
    output:
        json = "{tool}/{vers}/{dtbs}/{args}.build.size.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "database": wildcards.dtbs,
                                    "arguments": str2args(wildcards.args),
                                    "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers][wildcards.dtbs]["fixed_args"])}
    run:
        out_json = json_default(mode="build", category="size", config=params.config)
        out_json["metrics"]["size_bytes"] = {}
        out_json["metrics"]["size_bytes"]["total"] = 0
        out_json["metrics"]["size_bytes"]["files"] = {}
        with open(input.fsize, "r") as file:
            for line in file:
                s, f = line.rstrip().split("\t")
                out_json["metrics"]["size_bytes"]["total"] += int(s)
                out_json["metrics"]["size_bytes"]["files"][f] = int(s)
        json_write(out_json, output.json)
        
