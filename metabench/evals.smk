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
        evals_classify = expand("{i}.{ext}", i=input_evals_classify(), ext=["classify.stats.json", 
                                                                            "classify.evals.json",
                                                                            "profile.stats.json",
                                                                            "profile.evals.json"])


rule stats_profile:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.bioboxes"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.stats.json"
    params:
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool,
                                                           "version": wildcards.vers,
                                                           "sample": wildcards.samp,
                                                           "database": wildcards.dtbs,
                                                           "database_arguments": str2args(wildcards.dtbs_args),
                                                           "arguments": str2args(wildcards.args)}),
        ranks = " ".join(config["ranks"])
    shell: 
        """
        stats=$(grep -v "^@" {input.bioboxes} | awk 'BEGIN{{FS="\\t"}}{{sum_rank[$2]+=$5}}END{{for(r in sum_rank){{ print "    \\""r"\\": " sum_rank[r] ","}}}}')

echo "{{
{params.json_wildcards}\\"stats\\": {{
${{stats::-1}}
}}
}}" > {output.json}
        """

rule stats_classify:
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

rule evals_classify:
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
        gt = lambda wildcards: config["samples"][wildcards.samp]["classify"]
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals_classify.py \
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

rule evals_profile:
    input:
        bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.bioboxes"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.evals.json"
    log:
        "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.evals.log"
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
        gt = lambda wildcards: config["samples"][wildcards.samp]["profile"],
        threhsold_profile = " ".join(map(str,config["threhsold_profile"]))
        
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals_profile.py \
                --ranks {params.ranks} \
                --input-results {input.bioboxes} \
                --input-ground-truth {params.gt} \
                {params.db_profile} \
                --taxonomy {config[taxonomy]} \
                --taxonomy-files {params.taxonomy_files} \
                --output-json {output.json} \
                --thresholds {params.threhsold_profile} 2> {log}

echo "{{
{params.json_wildcards}\\"evals\\":
$(cat {output.json})
}}" > {output.json}
        """
