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
        srcdir("../envs/ganon_env.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db",
        db_folder = lambda wildcards: ("--input " + config["dbs"][wildcards.dtbs]["folder"]) if "folder" in config["dbs"][wildcards.dtbs] else "",
        db_extension = lambda wildcards: ("--input-extension " + config["dbs"][wildcards.dtbs]["extension"]) if "extension" in config["dbs"][wildcards.dtbs] else "",
        db_taxonomy = lambda wildcards: ("--taxonomy " + config["dbs"][wildcards.dtbs]["taxonomy"]) if "taxonomy" in config["dbs"][wildcards.dtbs] else "",
        db_taxonomy_files = lambda wildcards: ("--taxonomy-files " + config["dbs"][wildcards.dtbs]["taxonomy_files"]) if "taxonomy_files" in config["dbs"][wildcards.dtbs] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        {params.path}ganon build-custom \
                           --db-prefix {params.outprefix} \
                           {params.db_folder} \
                           {params.db_extension} \
                           {params.db_taxonomy} \
                           {params.db_taxonomy_files} \
                           --threads {threads} \
                           --verbose \
                           {params.fixed_args} \
                           {params.args} > {log} 2>&1

        if [[ -e "{params.outprefix}.hibf" ]]; then
            ln -sr "{params.outprefix}.hibf" "{params.outprefix}.ibf"
        fi
        """

rule ganon_build_size:
    input:
        db1 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.ibf",
        db2 = "ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db.tax"
    output:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "du --bytes --dereference {input} > {output}"  # output in bytes


rule ganon_binning:
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
        srcdir("../envs/ganon_env.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        outprefix = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}", 
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args),
        reassign = lambda wildcards: 1 if "--reassign" in wildcards.b_args else 0 
    shell:
        """
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}ganon classify \
                               --db-prefix {params.dbprefix} \
                               --single-reads {input.fq1} \
                               --output-prefix {params.outprefix} \
                               --threads {threads} \
                               --verbose --output-lca \
                               {params.fixed_args} \
                               {params.args} > {log} 2>&1
        else # paired-end
            {params.path}ganon classify \
                               --db-prefix {params.dbprefix} \
                               --paired-reads {input.fq1} {params.input_fq2} \
                               --output-prefix {params.outprefix} \
                               --threads {threads} \
                               --verbose --output-lca \
                               {params.fixed_args} \
                               {params.args} > {log} 2>&1
        fi

        if [[ "{params.reassign}" -eq "1" ]]; then
            mv {params.outprefix}.all {params.outprefix}.lca
        fi 
        """

rule ganon_binning_format:
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

        # Check with .tax entries at "assembly" or "sequence" level and report parent taxa
        # also report entries not matching .tax (-a 1)
        # Check if end of read id is "/1" and remove it

        join -1 2 -2 1 -t$'\t' -o "1.1,1.2,2.2,2.3" -a 1 \
        <(sort -t$'\t' -k 2,2 {input.lca}) \
        <(sort -t$'\t' -k 1,1 {input.dbtax}) | \
        awk 'BEGIN{{FS=OFS="\t"}}
            {{
            if(substr($1,length($1)-1)=="/1"){{
                $1=substr($1,0,length($1)-2);
            }};
            if($4=="assembly"){{
                print $1,$3,$2,"";
            }}else{{
                if($4=="sequence"){{
                    print $1,$3,"",$2;
                }}else{{
                    print $1,$2,"","";
                }}
            }}
            }}' >> {output.bioboxes}
        """


rule ganon_profiling:
    input:
        rep = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep"
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    conda:
        srcdir("../envs/ganon_env.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        {params.path}ganon report \
                           --db-prefix {params.dbprefix} \
                           --input {input.rep} \
                           --output-prefix {output.bioboxes} \
                           --output-format bioboxes {params.args} > {log} 2>&1
        mv {output.bioboxes}.tre {output.bioboxes}
        """