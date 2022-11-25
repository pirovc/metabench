import glob
import os
from collections import defaultdict

workdir: config["workdir"]
include: "util.py"

def input_all():
    out = []
    for binning_file in glob.glob('**/*.binning.bioboxes', recursive=True):
        # without ".bioboxes"
        binning_file_prefix = os.path.splitext(binning_file)[0]
        # Just generate for prefixes with sample in config file
        sample = binning_file_prefix.split("/")[2]
        if sample in config["samples"].keys():
            out.append(binning_file_prefix)

    for profiling_file in glob.glob('**/*.profiling.bioboxes', recursive=True):
        # without ".bioboxes"
        profiling_file_prefix = os.path.splitext(profiling_file)[0]
        # Just generate for prefixes with sample in config file
        sample = profiling_file_prefix.split("/")[2]
        if sample in config["samples"].keys():
            out.append(profiling_file_prefix)

    return out

rule all:
    input:
        expand("{i}.{ext}", i=input_all(), ext=["evals.json"])


rule evals_binning_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
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
        gt = lambda wildcards: config["samples"][wildcards.samp]["binning"]
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
                --output-cumulative {output.cumu_json} \
                --output-rank {output.rank_json} 2> {log}
        """

rule evals_binning:
    input:
        cumu_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.cumu.json",
        rank_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.rank.json"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.evals.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args)}
    run:
        out_json = json_default(report="binning", category="evals", config=params.config)
        out_json["metrics"] = {}
        out_json["metrics"]["cumulative-based"] = json_load(input.cumu_json)
        out_json["metrics"]["rank-based"] = json_load(input.rank_json)
        json_write(out_json, output.json)


rule evals_profiling_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
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
        threhsold_profiling = " ".join(map(str,config["threhsold_profiling"]))
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
        json_tmp = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.evals.tmp.json"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.evals.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "binning_arguments": str2args(wildcards.b_args),
                                    "profiling_arguments": str2args(wildcards.p_args)}
    run:
        out_json = json_default(report="profiling", category="evals", config=params.config)
        out_json["metrics"] = json_load(input.json_tmp)
        json_write(out_json, output.json)

