import glob
import os
from collections import defaultdict

workdir: config["workdir"]
include: "util.py"

def input_all(t):
    out = []
    suffix_len = len(".bioboxes.gz")
    for file in glob.glob('**/*.' + t + '.bioboxes.gz', recursive=True):
        # without ".bioboxes.gz"
        file_prefix = file[:-suffix_len]
        # Just generate for prefixes with sample in config file
        sample = file_prefix.split("/")[2]
        if sample in config["samples"].keys():
            out.append(file_prefix)
    return out

rule all:
    input:
        binning = expand("{i}.{ext}", i=input_all("binning"), ext=["bench.json", "updated_json"]),
        profiling = expand("{i}.{ext}", i=input_all("profiling"), ext=["bench.json", "updated_json"]),
        
rule evals_binning_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes.gz"
    output:
        cumu_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.cumu.json"),
        rank_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.rank.json"),
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.log"
    params:
        scripts_path = srcdir("../scripts/"),
        ranks = " ".join(config["ranks"]),
        taxonomy_files = config["taxonomy_files"],
        db_profile = lambda wildcards: "--input-database-profile " + config["dbs"][wildcards.dtbs] if "dbs" in config and wildcards.dtbs in config["dbs"] else "",
        gt = lambda wildcards: config["samples"][wildcards.samp]["binning"],
        threhsold_binning = " ".join(map(str,config["threhsold_binning"])) if "threhsold_binning" in config else "0"
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals_binning.py \
                --ranks {params.ranks} \
                --input-results {input.bioboxes} \
                --input-ground-truth {params.gt} \
                {params.db_profile} \
                --taxonomy {config[taxonomy]} \
                --taxonomy-files {params.taxonomy_files} \
                --output-rank {output.rank_json} \
                --output-cumu {output.cumu_json} \
                --thresholds {params.threhsold_binning} 2> {log}
        """

rule evals_binning:
    input:
        bench_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.json",
        cumu_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.cumu.json",
        rank_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.rank.json"
    output:
        flag_update = touch("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.updated_json")
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args)}
    run:
        out_json = json_load(input.bench_json)
        out_json["evals"] = {} # delete old evals in case of rewrite
        out_json["evals"].update(json_load(input.cumu_json))
        out_json["evals"].update(json_load(input.rank_json))
        json_write(out_json, input.bench_json)


rule evals_profiling_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes.gz"
    output:
        json_tmp = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.evals.tmp.json")
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.evals.log"
    params:
        scripts_path = srcdir("../scripts/"),
        ranks = " ".join(config["ranks"]),
        taxonomy_files = config["taxonomy_files"],
        db_profile = lambda wildcards: "--input-database-profile " + config["dbs"][wildcards.dtbs] if "dbs" in config and wildcards.dtbs in config["dbs"] else "",
        gt = lambda wildcards: config["samples"][wildcards.samp]["profiling"],
        threhsold_profiling = " ".join(map(str,config["threhsold_profiling"]))  if "threhsold_profiling" in config else "0"
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals_profiling.py \
                --ranks {params.ranks} \
                --input-results {input.bioboxes} \
                --input-ground-truth {params.gt} \
                {params.db_profile} \
                --taxonomy {config[taxonomy]} \
                --taxonomy-files {params.taxonomy_files} \
                --output-json {output.json_tmp} \
                --thresholds {params.threhsold_profiling} 2> {log}
        """

rule evals_profiling:
    input:
        bench_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.json",
        json_tmp = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.evals.tmp.json"
    output:
        flag_update = touch("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.updated_json")
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args),
                                    "profiling_arguments": str2args(wildcards.p_args)}
    run:
        out_json = json_load(input.bench_json)
        out_json["evals"] = {}  # delete old evals in case of rewrite
        out_json["evals"] = json_load(input.json_tmp)
        json_write(out_json, input.bench_json)

