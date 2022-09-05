import numpy as np
from itertools import product

workdir: config["workdir"]
include: "util.py"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_all(wildcards, ext: list):
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                for dtbs in config["dbs"]:
                    # define path to output
                    path = tool + "/" + vers + "/" + dtbs + "/"
                    # build product of all parameters (single or range)
                    for e in ext:
                        for p in params_product(config["run"][tool]["params"]):
                            out.append(path + join_params(p) + "." + e)
    
    return out


rule all:
    input: lambda wildcards: unpack(input_all(wildcards, ext=["build.bench.json", "size.json"]))

rule time:
    input: bench = "{tool}/{vers}/{dtbs}/{prms}.build.bench.tsv"
    output: json = "{tool}/{vers}/{dtbs}/{prms}.build.bench.json"
    params: 
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)}),
        json_bench = lambda wildcards, input: json_bench(input.bench),
    shell: 
        """
echo "{{
{params.json_wildcards}\\"bench\\": {{
{params.json_bench}
}}
}}" > {output.json}
        """

rule fsize:
    input: fsize = "{tool}/{vers}/{dtbs}/{prms}.size.tsv"
    output: json = "{tool}/{vers}/{dtbs}/{prms}.size.json"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)})
    shell: 
        """
        # index size bytes
        file_size_bytes=$(awk -F "\\t" '{{print "    \\"" $2 "\\": " $1 ","}}' {input.fsize})
        total_size_bytes=$(awk -F "\\t" '{{s+=$1}}END{{print s}}' {input.fsize})

echo "{{
{params.json_wildcards}\\"size\\":
{{
${{file_size_bytes}}
    \\"total_size_bytes\\": ${{total_size_bytes}}
}}
}}" > {output.json}
        """