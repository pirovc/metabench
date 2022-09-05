import numpy as np
from itertools import product

workdir: config["workdir"]
include: "util.py"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_all(wildcards, ext: list):
    out = []
    for tool in config["tools"]:
        for vers in config["tools"][tool]:
            for samp in config["samples"]:
                if tool in config["run"]:
                    for dtbs in config["run"][tool]["dbs"]:
                        # define path to output
                        path = tool + "/" + vers + "/" + samp + "/" + dtbs + "/"
                        # build product of all parameters (single or range)
                        pparams = params_product(config["run"][tool]["params"])
                        for e in ext:
                            for p in pparams:
                                out.append(path + join_params(p) + "." + e)

    return out


rule all:
    input: lambda wildcards: unpack(input_all(wildcards, ext=["bioboxes", "classify.bench.json"]))


rule time:
    input: bench = "{tool}/{vers}/{samp}/{dtbs}/{prms}.classify.bench.tsv"
    output: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.classify.bench.json"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "sample": wildcards.samp, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)}),
            json_bench = lambda wildcards, input: json_bench(input.bench),
    shell: 
        """
echo "{{
{params.json_wildcards}\\"bench\\": {{
{params.json_bench}
}}
}}" > {output.json}
        """