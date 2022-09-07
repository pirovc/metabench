workdir: config["workdir"]
include: "util.py"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_classify():
    out = []
    for tool in config["tools"]:
        if tool in config["run"]:
            for vers in config["tools"][tool]:
                if vers in config["run"][tool]:
                    for samp in config["samples"]:
                        for dtbs in config["run"][tool][vers]["dbs"]:
                            # Get all folders in database prefix
                            db_dir = config["run"][tool][vers]["dbs"][dtbs]
                            for dtbs_args in [d for d in os.listdir(db_dir) if os.path.isdir(os.path.join(db_dir, d))]:
                                # define path to output
                                path = tool + "/" + vers + "/" + samp + "/" + dtbs + "/" + dtbs_args + "/"
                                # build product of all arguments (single or range)
                                for args in args_product(config["run"][tool][vers]["args"]):
                                    out.append(path + join_args(args))
    return out

def input_profile():
    out = []
    for prefix in input_classify():
        tool = prefix.split("/")[0]
        vers = prefix.split("/")[1]
        if "args_profile" in config["run"][tool][vers]:
            for args in args_product(config["run"][tool][vers]["args_profile"]):
                out.append(prefix + "/" + join_args(args))
    return out

rule all:
    input:
        classify = expand("{i}.{ext}", i=input_classify(), ext=["classify.bioboxes", "classify.bench.json"]),
        profile = expand("{i}.{ext}", i=input_profile(), ext=["profile.bioboxes"])

rule time:
    input:
        bench = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.bench.tsv"
    output:
        json = "{tool}/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.classify.bench.json"
    params:
        json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, 
                                                           "version": wildcards.vers,
                                                           "sample": wildcards.samp,
                                                           "database": wildcards.dtbs,
                                                           "database_arguments": str2args(wildcards.dtbs_args),
                                                           "arguments": str2args(wildcards.args),
                                                           "fixed_arguments": dict2args(config["run"][wildcards.tool][wildcards.vers]["fixed_args"])}),
        json_bench = lambda wildcards, input: json_bench(input.bench),
    shell: 
        """
echo "{{
{params.json_wildcards}\\"bench\\": {{
{params.json_bench}
}}
}}" > {output.json}
        """