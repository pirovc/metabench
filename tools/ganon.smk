rule ganon_build:
    output: db1 = "ganon/{vers}/{dtbs}/{prms}.ibf",
            db2 = "ganon/{vers}/{dtbs}/{prms}.tax"
    benchmark: "ganon/{vers}/{dtbs}/{prms}.build.bench.tsv"
    log: "ganon/{vers}/{dtbs}/{prms}.log"
    threads: config["threads"]
    conda: srcdir("../envs/ganon.yaml")
    params: folder = lambda wildcards: os.path.abspath(config["dbs"][wildcards.dtbs]["folder"]),
            extens = lambda wildcards: config["dbs"][wildcards.dtbs]["extension"],
            outprefix = "ganon/{vers}/{dtbs}/{prms}",
            path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
            fixed_params = lambda wildcards: config["run"]["ganon"]["fixed_params"] if "fixed_params" in config["run"]["ganon"] else "",
            params = lambda wildcards: str2params(wildcards.prms),
            taxonomy_files = " ".join(config["taxonomy_files"]) if "taxonomy_files" in config else ""
    shell: 
        """
        # if path is provided, deactivate conda
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        {params.path}ganon build-custom \
        --db-prefix {params.outprefix} \
        --input {params.folder} \
        --input-extension {params.extens} \
        --threads {threads} \
        --taxonomy {config[taxonomy]} \
        --taxonomy-files {params.taxonomy_files} \
        {params.fixed_params} {params.params} > {log} 2>&1
        """

rule ganon_build_size:
    input: db1 = "ganon/{vers}/{dtbs}/{prms}.ibf",
           db2 = "ganon/{vers}/{dtbs}/{prms}.tax"
    output: "ganon/{vers}/{dtbs}/{prms}.size.tsv"
    shell: "du --block-size=1 {input} > {output}" #output in bytes


rule ganon_classify:
    input: fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output: rep=temp("ganon/{vers}/{samp}/{dtbs}/{prms}.rep"),
            lca=temp("ganon/{vers}/{samp}/{dtbs}/{prms}.lca"),
            tre=temp("ganon/{vers}/{samp}/{dtbs}/{prms}.tre")
    benchmark: "ganon/{vers}/{samp}/{dtbs}/{prms}.classify.bench.tsv"
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
            {params.path}ganon classify --verbose --output-lca --db-prefix {params.dbprefix} -s {input.fq1} -t {threads} {params.fixed_params} {params.params} -o {params.outprefix} > {log} 2>&1
        else # paired-end
            {params.path}ganon classify --verbose --output-lca--db-prefix {params.dbprefix} -p {input.fq1} {params.input_fq2} -t {threads} {params.fixed_params} {params.params} -o {params.outprefix} > {log} 2>&1
        fi
        """

rule ganon_classify_format:
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
