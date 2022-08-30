import glob
import os

workdir: config["workdir"]
include: "util.py"


rule all:
    input: stats = expand("{prefix}.stats.json", prefix=[os.path.splitext(file)[0] for file in glob.glob('**/*.bioboxes', recursive=True)])


rule stats:
    input: bioboxes = "{tool}/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    output: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.stats.json"
    log: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.stats.log"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "sample": wildcards.samp, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)}),
            scripts_path = srcdir("../scripts/"),
            ranks = " ".join(config["ranks"]),
            taxonomy_files = " ".join(config["taxonomy_files"])
    conda: srcdir("../envs/eval.yaml")
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
