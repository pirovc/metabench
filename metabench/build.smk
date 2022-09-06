workdir: config["workdir"]
include: "util.py"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_all(wildcards, ext: list):
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                if vers in config["run"][tool]:
                    for dtbs in config["dbs"]:
                        # define path to output
                        path = tool + "/" + vers + "/" + dtbs + "/"
                        # build product of all arguments (single or range)
                        for args in args_product(config["run"][tool][vers][dtbs]["args"]):
                            # For every final final extension
                            for e in ext:
                                out.append(path + join_args(args) + "." + e)
    
    return out


rule all:
    input:
        lambda wildcards: unpack(input_all(wildcards, ext=["build.bench.json", "size.json"]))

rule time:
    input:
        bench = "{tool}/{vers}/{dtbs}/{args}.build.bench.tsv"
    output:
        json = "{tool}/{vers}/{dtbs}/{args}.build.bench.json"
    params: 
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "database": wildcards.dtbs, "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers][wildcards.dtbs]["fixed_args"]), "arguments": str2args(wildcards.args)}),
        json_bench = lambda wildcards, input: json_bench(input.bench)
    shell: 
        """
echo "{{
{params.json_wildcards}\\"bench\\": {{
{params.json_bench}
}}
}}" > {output.json}
        """

rule fsize:
    input:
        fsize = "{tool}/{vers}/{dtbs}/{args}.size.tsv"
    output:
        json = "{tool}/{vers}/{dtbs}/{args}.size.json"
    params:
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "database": wildcards.dtbs, "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers][wildcards.dtbs]["fixed_args"]), "arguments": str2args(wildcards.args)})
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