rule ganon_build:
    output:
        db1 = "ganon/{vers}/{dtbs}/{args}/ganon_db.ibf",
        db2 = "ganon/{vers}/{dtbs}/{args}/ganon_db.tax"
    benchmark:
        repeat("ganon/{vers}/{dtbs}/{args}.build.bench.tsv", config["repeat"])
    log:
        "ganon/{vers}/{dtbs}/{args}.build.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/ganon.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{dtbs}/{args}/ganon_db",
        db = lambda wildcards: config["dbs"][wildcards.dtbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.args)
    shell: 
        """
        # if path is provided, deactivate conda
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        {params.path}ganon build-custom \
        --db-prefix {params.outprefix} \
        --input {params.db[folder]} \
        --input-extension {params.db[extension]} \
        --taxonomy {params.db[taxonomy]} \
        --taxonomy-files {params.db[taxonomy_files]} \
        --threads {threads} \
        --verbose \
        {params.fixed_args} {params.args} > {log} 2>&1
        """

rule ganon_build_size:
    input:
        db1 = "ganon/{vers}/{dtbs}/{args}/ganon_db.ibf",
        db2 = "ganon/{vers}/{dtbs}/{args}/ganon_db.tax"
    output:
        "ganon/{vers}/{dtbs}/{args}.build.size.tsv"
    shell:
        "du --block-size=1 {input} > {output}"  # output in bytes


rule ganon_classify:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output:
        lca=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.lca"),
        rep=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.rep"),
        tre=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.tre")
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.bench.tsv", config["repeat"])
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/ganon.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}", 
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.args)
    shell:
        """
        # if path is provided, deactivate conda
        # if [[ ! -z "{params.path}" ]]; then
        #     source deactivate;
        # fi
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}ganon classify \
                               --db-prefix {params.dbprefix} \
                               --single-reads {input.fq1} \
                               --output-prefix {params.outprefix} \
                               --threads {threads} \
                               --verbose --output-lca \
                               {params.fixed_args} {params.args} > {log} 2>&1
        else # paired-end
            {params.path}ganon classify \
                               --db-prefix {params.dbprefix} \
                               --paired-reads {input.fq1} {params.input_fq2} \
                               --output-prefix {params.outprefix} \
                               --threads {threads} \
                               --verbose --output-lca \
                               {params.fixed_args} {params.args} > {log} 2>&1
        fi
        """

rule ganon_classify_format:
    input: 
        lca="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.lca",
        dbtax = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.binning.bioboxes"
    params:
        header = lambda wildcards: header_bioboxes_binning("ganon", wildcards)
    shell:
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # if numeric (taxid) or assembly
        awk 'FS="\\t"{{if($2 ~ /^[0-9]+$/){{print substr($1,1,length($1)-2)"\\t0\\t"$2}}}}' {input.lca} >> {output.bioboxes}

        # get taxid if classified at assembly level
        join -1 2 -2 1 \
        <(awk 'FS="\\t"{{if( $2 !~ /^[0-9]+$/){{print substr($1,1,length($1)-2)"\\t"$2"\\t"$3}}}}' {input.lca} | sort -k 2,2) \
        <(cut -f 1,2 {input.dbtax} | sort | uniq | sort -k 1,1) \
        -t$'\\t' -o "1.1,1.2,2.2" >> {output.bioboxes}
        """

rule ganon_profiling_format:
    input:
        tre = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.tre"
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profiling.bioboxes.log"
    conda:
        srcdir("../envs/evals.yaml")
    params:
        scripts_path = srcdir("../scripts/"),
        ranks = lambda wildcards: " ".join(config["ranks_profiling"]),
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"],
        header = lambda wildcards: header_bioboxes_profiling("ganon",
                                                           config["ranks_profiling"],
                                                           [os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"],
                                                           wildcards),
    shell: 
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}
        python3 {params.scripts_path}profile.py \
                                    --taxid-col 1 \
                                    --perc-col 8 \
                                    --input-file {input.tre} \
                                    --taxonomy custom \
                                    --taxonomy-files {params.taxonomy_files} \
                                    --ranks {params.ranks} >> {output.bioboxes} 2> {log}

        """
