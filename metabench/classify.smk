import numpy as np
from itertools import product

workdir: config["workdir"]
include: "util.py"
include: "../tools/ganon.smk"
include: "../tools/kraken2.smk"


def input_all(wildcards, ext: list):
    out = []
    for tool in config["tools"]:
        for vers in config["tools"][tool]:
            for samp in config["samples"]:
                if tool in config["run"]:
                    for run in config["run"][tool]:
                        for dtbs in config["run"][tool]["dbs"]:
                            # define path to output
                            path = tool + "/" + vers + "/" + samp + "/" + dtbs + "/"

                            # define filename based on product of parameters
                            params = []
                            # ranged
                            if "params" in config["run"][tool]:
                                for p, v in config["run"][tool]["params"].items():
                                    if isinstance(v, list):
                                        rparams = []
                                        for range_val in list(np.arange(v[0], v[1], v[2])) + [v[1]]: # include stop
                                            rparams.append(params2str(p, round(range_val, 2)))  # round to 2 decimal
                                        params.append(tuple(rparams))
                                    else:
                                        params.append(tuple([params2str(p,v)]))

                            for e in ext:
                                if params:
                                    for p in product(*params):
                                        out.append(path + join_params(p) + "." + e)
                                else:
                                    out.append(path + "default." + e)
    return out


rule all:
    input: lambda wildcards: unpack(input_all(wildcards, ext=["bioboxes", "time.json"]))


rule time:
    input: time = "{tool}/{vers}/{samp}/{dtbs}/{prms}.time"
    output: json = "{tool}/{vers}/{samp}/{dtbs}/{prms}.time.json"
    params: json_wildcards = lambda wildcards: json_wildcards({"tool": wildcards.tool, "version": wildcards.vers, "sample": wildcards.samp, "database": wildcards.dtbs, "parameters": str2params(wildcards.prms)})
    shell: 
        """
        # elapsed time
        timeelapsed=$(grep -oP "(?<=Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ).+" {input.time})
        # time in seconds
        timesec=$(echo ${{timeelapsed}} | awk '{{l=split($0,a,":");if(l==2){{sec=(a[1]*60)+a[2]}}else{{sec=(a[1]*3600)+(a[2]*60)+a[3]}};print sec }}')
        # peak memory in kilobytes 
        memkbytes=$(grep -oP "(?<=Maximum resident set size \\(kbytes\\): ).+" {input.time})
        # peak memory in bytes
        membytes=$((memkbytes*1000))

echo "{{
{params.json_wildcards}\\"time\\": {{
    \\"time_seconds\\": ${{timesec}},
    \\"mem_bytes\\": ${{membytes}}
}}
}}" > {output.json}
        """