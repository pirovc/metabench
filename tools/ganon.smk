rule ganon_classify:
    input: fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
           dbibf = lambda wildcards: os.path.abspath(config["run"]["ganon"]["dbs"][wildcards.dtbs] + ".ibf"),
           dbtax = lambda wildcards: os.path.abspath(config["run"]["ganon"]["dbs"][wildcards.dtbs] + ".tax")
    output: rep="ganon/{vers}/{samp}/{dtbs}/{prms}.rep",
            lca="ganon/{vers}/{samp}/{dtbs}/{prms}.lca",
            time="ganon/{vers}/{samp}/{dtbs}/{prms}.time",
    log: "ganon/{vers}/{samp}/{dtbs}/{prms}.log"
    threads: config["threads"]
    conda: srcdir("../envs/ganon.yaml")
    params: outprefix = "ganon/{vers}/{samp}/{dtbs}/{prms}", 
            path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
            dbprefix = lambda wildcards: config["run"]["ganon"]["dbs"][wildcards.dtbs],
            input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
            fixed_params = lambda wildcards: config["run"]["ganon"]["fixed_params"] if "fixed_params" in config["run"]["ganon"] else "",
            params = lambda wildcards: str2params(wildcards.prms)
    shell:
        """
        # if path is provided, deactivate conda
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            /usr/bin/time -v --output={output.time} {params.path}ganon classify --verbose --output-lca --db-prefix {params.dbprefix} -s {input.fq1} -t {threads} {params.fixed_params} {params.params} -o {params.outprefix} > {log} 2>&1
        else # paired-end
            /usr/bin/time -v --output={output.time} {params.path}ganon classify --verbose --output-lca--db-prefix {params.dbprefix} -p {input.fq1} {params.input_fq2} -t {threads} {params.fixed_params} {params.params} -o {params.outprefix} > {log} 2>&1
        fi
        """

rule ganon_format:
    input: lca="ganon/{vers}/{samp}/{dtbs}/{prms}.lca",
           dbtax = lambda wildcards: os.path.abspath(config["run"]["ganon"]["dbs"][wildcards.dtbs] + ".tax")
    output: bioboxes = "ganon/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    shell:
        """
        # bioboxes header
        printf "@Version:0.9.1\\n@SampleID:ganon {wildcards.vers} {wildcards.samp} {wildcards.dtbs} {wildcards.prms}\\n@@SEQUENCEID\\tBINID\\tTAXID\\n" > {output.bioboxes}

        # if numeric (taxid) or assembly
        awk 'FS="\\t"{{if($2 ~ /^[0-9]+$/){{print substr($1,1,length($1)-2)"\\t0\\t"$2}}}}' {input.lca} >> {output.bioboxes}

        # get taxid if classified at assembly level
        join -1 2 -2 1 \
        <(awk 'FS="\\t"{{if( $2 !~ /^[0-9]+$/){{print substr($1,1,length($1)-2)"\\t"$2"\\t"$3}}}}' {input.lca} | sort -k 2,2) \
        <(cut -f 1,2 {input.dbtax} | sort | uniq | sort -k 1,1) \
        -t$'\\t' -o "1.1,1.2,2.2" >> {output.bioboxes}
        """
