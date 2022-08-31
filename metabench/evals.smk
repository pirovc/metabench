import glob
import os

workdir: config["workdir"]
include: "util.py"

def get_results_prefix():
    return [os.path.splitext(file)[0] for file in glob.glob('**/*.bioboxes', recursive=True)]

def get_results_prefix_filtered():
    res = get_results_prefix()
    filtered_res = []
    for samp in config["samples"].keys():
        for dtbs in config["dbs"].keys():
            for r in res:
                if samp in r and dtbs in r:
                    filtered_res.append(r)
    return filtered_res


rule all:
    input: stats = expand("{prefix}.stats.json", prefix=get_results_prefix()),
           evals = expand("{prefix}.evals.json", prefix=get_results_prefix_filtered())


rule stats:
    input: bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    output: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.stats.json"
    log: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.stats.log"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "sample": wildcards.samp, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)}),
            scripts_path = srcdir("../scripts/"),
            ranks = " ".join(config["ranks"]),
            taxonomy_files = " ".join(config["taxonomy_files"])
    conda: srcdir("../envs/evals.yaml")
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
    input: bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    output: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.evals.json",
            cumu_json = temp("{tool}/{vers}/{samp}/{dtbs}/{prms}.evals.cumu.json"),
            rank_json = temp("{tool}/{vers}/{samp}/{dtbs}/{prms}.evals.rank.json"),
    log: "{tool}/{vers}/{samp}/{dtbs}/{prms}.evals.log"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "sample": wildcards.samp, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)}),
            scripts_path = srcdir("../scripts/"),
            ranks = " ".join(config["ranks"]),
            taxonomy_files = " ".join(config["taxonomy_files"]),
            db_profile = lambda wildcards: config["dbs"][wildcards.dtbs],
            gt = lambda wildcards: config["samples"][wildcards.samp]
    conda: srcdir("../envs/evals.yaml")
    shell: 
        """
        python3 {params.scripts_path}evals.py \
                --ranks {params.ranks} \
                --input-results {input.bioboxes} \
                --input-ground-truth {params.gt} \
                --input-database-profile {params.db_profile} \
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
