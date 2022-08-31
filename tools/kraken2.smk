rule kraken2_classify:
    input: fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
           dbfolder = lambda wildcards: os.path.abspath(config["run"]["kraken2"]["dbs"][wildcards.dtbs])
    output: res=temp("kraken2/{vers}/{samp}/{dtbs}/{prms}.res"),
            time="kraken2/{vers}/{samp}/{dtbs}/{prms}.time"
    log: "kraken2/{vers}/{samp}/{dtbs}/{prms}.log"
    threads: config["threads"]
    conda: srcdir("../envs/kraken2.yaml")
    params: outprefix = "kraken2/{vers}/{samp}/{dtbs}/{prms}",
            path = lambda wildcards: config["tools"]["kraken2"][wildcards.vers],
            input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
            fixed_params = lambda wildcards: config["run"]["kraken2"]["fixed_params"] if "fixed_params" in config["run"]["kraken2"] else "",
            params = lambda wildcards: str2params(wildcards.prms)
    shell:
        """
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            /usr/bin/time -v --output={output.time} {params.path}kraken2 --db {input.dbfolder} --threads {threads} --output {output.res} --gzip-compressed {params.params} {input.fq1} > {log} 2>&1
        else # paired-end
            /usr/bin/time -v --output={output.time} {params.path}kraken2 --db {input.dbfolder} --threads {threads} --output {output.res} --gzip-compressed {params.params} --paired {input.fq1} {params.input_fq2} > {log} 2>&1
        fi
        """

rule kraken2_format:
    input: res="kraken2/{vers}/{samp}/{dtbs}/{prms}.res",
    output: "kraken2/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    params: input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
    shell:
        """
        # bioboxes header
        printf "@Version:0.9.1\\n@SampleID:kraken2 {wildcards.vers} {wildcards.samp} {wildcards.dtbs} {wildcards.prms}\\n@@SEQUENCEID\\tBINID\\tTAXID\\n" > {output}

        # output header changes when sinle or paired (/1 or nothing)
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            header_suffix=2;
        else # paired-end
            header_suffix=0;
        fi

        grep "^C" {input.res} | awk -v header_suffix="${{header_suffix}}" 'FS="\\t"{{print substr($2,1,length($2)-header_suffix)"\\t0\\t"$3}}' >> {output}
        """
