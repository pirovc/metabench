import glob
import os

workdir: config["workdir"]
include: "util.py"

def input_evals_classify():
    out = []
    for file in glob.glob('**/*.classify.bioboxes', recursive=True):
        # without ".classify.bioboxes"
        prefix = os.path.splitext(os.path.splitext(file)[0])[0]
        # Just generate for prefixes with sample in config fige
        sample = prefix.split("/")[3]
        for sample in config["samples"].keys():
            out.append(prefix)
    return out

rule all:
    input:
        evals_classify = expand("{i}.{ext}", i=input_evals_classify(), ext=["classify.stats.json", "classify.evals.json"])


rule stats:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.bioboxes"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.stats.json"
    log:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.stats.log"
    params:
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool,
                                                           "version": wildcards.vers,
                                                           "sample": wildcards.samp,
                                                           "database": wildcards.dtbs,
                                                           "database_arguments": str2args(wildcards.dtbs_args),
                                                           "arguments": str2args(wildcards.args)}),
        scripts_path = srcdir("../scripts/"),
        ranks = " ".join(config["ranks"]),
        taxonomy_files = " ".join(config["taxonomy_files"])
    conda:
        srcdir("../envs/evals.yaml")
    shell: 
        """
        stats=$(python3 {params.scripts_path}stats.py \
                        --ranks {params.ranks} \
                        --input-file {input.bioboxes} \
                        --taxonomy {config[taxonomy]} \
                        --taxonomy-files {params.taxonomy_files} 2> {log})
echo "{{
{params.json_wildcards}\\"stats\\": ${{stats}}
}}" > {output.json}
        """

rule evals:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.bioboxes"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.evals.json",
        cumu_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.evals.cumu.json"),
        rank_json = temp("{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.evals.rank.json"),
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.evals.log"
    params:
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool,
                                                           "version": wildcards.vers,
                                                           "sample": wildcards.samp,
                                                           "database": wildcards.dtbs,
                                                           "database_arguments": str2args(wildcards.dtbs_args),
                                                           "arguments": str2args(wildcards.args)}),
        scripts_path = srcdir("../scripts/"),
        ranks = " ".join(config["ranks"]),
        taxonomy_files = " ".join(config["taxonomy_files"]),
        db_profile = lambda wildcards: "--input-database-profile " + config["dbs"][wildcards.dtbs] if "dbs" in config and wildcards.dtbs in config["dbs"] else "",
        gt = lambda wildcards: config["samples"][wildcards.samp]
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals.py \
                --ranks {params.ranks} \
                --input-results {input.bioboxes} \
                --input-ground-truth {params.gt} \
                {params.db_profile} \
                --taxonomy {config[taxonomy]} \
                --taxonomy-files {params.taxonomy_files} \
                --output-cumulative {output.cumu_json} \
                --output-rank {output.rank_json} 2> {log}

echo "{{
{params.json_wildcards}\\"evals\\":
{{
\\"cumulative\\":
$(cat {output.cumu_json})
,
\\"rank\\":
$(cat {output.rank_json})
}}
}}" > {output.json}
        """
