import numpy as np
from itertools import product

workdir: config["workdir"]
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
                                            rparams.append(p + "=" + str(round(range_val, 2)))  # round to 2 decimal
                                        params.append(tuple(rparams))
                                    else:
                                        params.append(tuple([p + "=" + str(v)]))

                            for e in ext:
                                if params:
                                    for p in product(*params):
                                        out.append(path + "_".join(p) + "." + e)
                                else:
                                    out.append(path + "default." + e)

    return out

rule all:
    input: lambda wildcards: unpack(input_all(wildcards, ext=["bioboxes", "stats"]))


rule stats:
    input: time = "{tool}/{vers}/{samp}/{dtbs}/{prms}.time"
    output: stats = "{tool}/{vers}/{samp}/{dtbs}/{prms}.stats"
    shell: 
        """
        #######
        # .stats: tool <tab> vers <tab> samp <tab> dtbs <tab> prms <tab> time elapsed <tab> time elapsed seconds <tab> peak memory bytes
        ######

        # elapsed time
        timeelapsed=$(grep -oP "(?<=Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ).+" {input.time})
        # time in seconds
        timesec=$(echo ${{timeelapsed}} | awk '{{l=split($0,a,":");if(l==2){{sec=(a[1]*60)+a[2]}}else{{sec=(a[1]*3600)+(a[2]*60)+a[3]}};print sec }}')

        # peak memory in kilobytes 
        memkbytes=$(grep -oP "(?<=Maximum resident set size \\(kbytes\\): ).+" {input.time})
        # peak memory in bytes
        membytes=$((memkbytes*1000))

        echo "{wildcards.tool}\t{wildcards.vers}\t{wildcards.samp}\t{wildcards.dtbs}\t{wildcards.prms}\t${{timeelapsed}}\t${{timesec}}\t${{membytes}}" > {output.stats}
        """