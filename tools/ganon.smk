rule ganon_build:
    output:
        db1 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.ibf",
        db2 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.tax"
    benchmark:
        repeat("ganon/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    log:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.log"
    threads:
        config["threads"]
    conda:
        #srcdir("../envs/ganon.yaml")
        srcdir("../envs/evals.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db",
        db = lambda wildcards: config["dbs"][wildcards.dtbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        # if path is provided, deactivate conda
        # if [[ ! -z "{params.path}" ]]; then
        #     source deactivate;
        # fi
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
        db1 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.ibf",
        db2 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.tax"
    output:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "du --block-size=1 {input} > {output}"  # output in bytes


rule ganon_classify:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output:
        lca=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.lca"),
        rep="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
        tre=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tre")
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    threads:
        config["threads"]
    conda:
        #srcdir("../envs/ganon.yaml")
        srcdir("../envs/evals.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}", 
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args)
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
        lca="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.lca",
        dbtax = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
    params:
        header = lambda wildcards: header_bioboxes_binning("ganon", wildcards)
    shell:
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # Check with .tax entries at "assembly" level and taxonomic level
        # also report entries not matching .tax (-a 1)
        # Check if end of read id is "/1" and remove it

        join -1 2 -2 1 <(sort -t$'\t' -k 2,2 {input.lca}) <(sort -t$'\t' -k 1,1 {input.dbtax}) -t$'\t' -o "1.1,1.2,2.2,2.3" -a 1 | awk 'BEGIN{{FS=OFS="\t"}}{{if(substr($1,length($1)-1)=="/1"){{$1=substr($1,0,length($1)-2)}};if($4=="assembly"){{print $1,$3,$2,""}}else{{if($4=="sequence"){{print $1,$3,"",$2}}else{{print $1,$3,"",""}}}}}}' >> {output.bioboxes}
        """


rule ganon_profiling:
    input:
        rep = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep"
    output:
        tre = temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.tre")
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    conda:
        #srcdir("../envs/ganon.yaml")
        srcdir("../envs/evals.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}", 
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        # if path is provided, deactivate conda
        # if [[ ! -z "{params.path}" ]]; then
        #     source deactivate;
        # fi
        
        {params.path}ganon report \
                           --db-prefix {params.dbprefix} \
                           --input {input.rep} \
                           --output-prefix {params.outprefix} {params.args} > {log} 2>&1
        """

rule ganon_profiling_format:
    input:
        tre = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.tre"
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes.log"
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
