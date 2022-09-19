import glob
import os
from collections import defaultdict

workdir: config["workdir"]
include: "util.py"

def input_all():
    out = []
    for file in glob.glob('**/*.binning.bioboxes', recursive=True):
        # without ".binning.bioboxes"
        prefix = os.path.splitext(os.path.splitext(file)[0])[0]
        # Just generate for prefixes with sample in config file
        sample = prefix.split("/")[2]
        if sample in config["samples"].keys():
            out.append(prefix)
    return out

rule all:
    input:
        expand("{i}.{ext}", i=input_all(), ext=["binning.stats.json", 
                                                "binning.evals.json",
                                                "profiling.stats.json",
                                                "profiling.evals.json"])

rule stats_profiling:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.bioboxes"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.stats.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "arguments": str2args(wildcards.args)}
    run:
        out_json = json_default(mode="profiling", category="stats", config=params.config)
        out_json["metrics"]["total_classified"] = {}
        out_json["metrics"]["total_classified"]["ranks"] = defaultdict(float)
        with open(input.bioboxes, "r") as file:
            for line in file:
                if line[0] == "@":
                    continue
                fields = line.rstrip().split("\t")
                # sum percentage of classification for each rank
                out_json["metrics"]["total_classified"]["ranks"][fields[1]] += float(fields[4])
        json_write(out_json, output.json)

rule stats_binning_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.bioboxes"
    output:
        json_tmp = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.stats.tmp.json")
    log:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.stats.log"
    conda:
        srcdir("../envs/evals.yaml")
    params:
         scripts_path = srcdir("../scripts/"),
         ranks = " ".join(config["ranks"]),
         taxonomy_files = config["taxonomy_files"]
    shell: 
        """
        python3 {params.scripts_path}stats_binning.py \
                --ranks {params.ranks} \
                --input-file {input.bioboxes} \
                --output-file {output.json_tmp} \
                --taxonomy {config[taxonomy]} \
                --taxonomy-files {params.taxonomy_files} 2> {log}
        """

rule stats_binning:
    input:
        json_tmp = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.stats.tmp.json"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.stats.json",
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "arguments": str2args(wildcards.args)}
    run:
        out_json = json_default(mode="binning", category="stats", config=params.config)
        out_json["metrics"] = json_load(input.json_tmp)
        json_write(out_json, output.json)


rule evals_binning_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.bioboxes"
    output:
        cumu_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.cumu.json"),
        rank_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.rank.json"),
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.log"
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
        cumu_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.cumu.json",
        rank_json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.rank.json"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.evals.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "arguments": str2args(wildcards.args)}
    run:
        out_json = json_default(mode="binning", category="evals", config=params.config)
        out_json["metrics"] = {}
        out_json["metrics"]["cumulative-based"] = json_load(input.cumu_json)
        out_json["metrics"]["rank-based"] = json_load(input.rank_json)
        json_write(out_json, output.json)


rule evals_profiling_script:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.bioboxes"
    output:
        json_tmp = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.evals.tmp.json")
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.evals.log"
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
        json_tmp = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.evals.tmp.json"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.evals.json"
    params:
        config = lambda wildcards: {"tool": wildcards.tool,
                                    "version": wildcards.vers,
                                    "sample": wildcards.samp,
                                    "database": wildcards.dtbs,
                                    "database_arguments": str2args(wildcards.dtbs_args),
                                    "arguments": str2args(wildcards.args)}
    run:
        out_json = json_default(mode="profiling", category="evals", config=params.config)
        out_json["metrics"] = json_load(input.json_tmp)
        json_write(out_json, output.json)

