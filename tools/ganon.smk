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
        srcdir("../envs/ganon-{vers}.yaml")
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
            rm -rf "{params.outprefix}.ibf" # remove if already existing
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
        "du --bytes --dereference --max-depth 0 {input} > {output}"  # output in bytes


rule ganon_binning:
    input:
        paired_reads = lambda wildcards: [os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
                                          os.path.abspath(config["samples"][wildcards.samp]["fq2"])],
        db = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.hibf"
    output:
        rep="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
        one=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one"),
        tre=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tre")
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    params:
        extra=lambda wildcards: str2args(wildcards.b_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
    wrapper:
        "file:/home/pirov/code/metabench/wrappers/ganon/{vers}/classify"


# rule ganon_binning2:
#     input:
#         fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
#     output:
#         one=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one"),
#         rep="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
#         tre=temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tre")
#     benchmark:
#         repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
#     log:
#         "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
#     threads:
#         config["threads"]
#     conda:
#         srcdir("../envs/ganon-{vers}.yaml")
#     params:
#         path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
#         outprefix = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}", 
#         dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
#         input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if "fq2" in config["samples"][wildcards.samp] else "",
#         fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"]),
#         args = lambda wildcards: str2args(wildcards.b_args)
#     shell:
#         """
#         if [[ -z "{params.input_fq2}" ]]; then # single-end
#             {params.path}ganon classify \
#                                --db-prefix {params.dbprefix} \
#                                --single-reads {input.fq1} \
#                                --output-prefix {params.outprefix} \
#                                --threads {threads} \
#                                --output-one \
#                                --verbose \
#                                {params.fixed_args} \
#                                {params.args} > {log} 2>&1
#         else # paired-end
#             {params.path}ganon classify \
#                                --db-prefix {params.dbprefix} \
#                                --paired-reads {input.fq1} {params.input_fq2} \
#                                --output-prefix {params.outprefix} \
#                                --threads {threads} \
#                                --output-one \
#                                --verbose \
#                                {params.fixed_args} \
#                                {params.args} > {log} 2>&1
#         fi

#         """

rule ganon_binning_format:
    input: 
        one="ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one",
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

        # 1.1 readid
        # 1.2 taxid
        # 2.2 parent taxid
        # 2.3 rank  

        join -1 2 -2 1 -t$'\t' -o "1.1,1.2,2.2,2.3" -a 1 \
        <(sort -t$'\t' -k 2,2 {input.one}) \
        <(sort -t$'\t' -k 1,1 {input.dbtax}) | \
        awk 'BEGIN{{FS=OFS="\t"}}
            {{
                if(substr($1,length($1)-1)=="/1"){{
                    $1=substr($1,0,length($1)-2);
                }};
                if($4=="assembly"){{
                    print $1,$3,$2,"";
                }}else if($4=="sequence"){{
                    print $1,$3,"",$2;
                }}else if($4=="file"){{
                    print $1,$3,"","";
                }}else{{
                    print $1,$2,"","";
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
        srcdir("../envs/ganon-{vers}.yaml")
    params:
        path = lambda wildcards: config["tools"]["ganon"][wildcards.vers],
        ranks = lambda wildcards: " ".join(config["default_ranks"]),
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db",
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        {params.path}ganon report \
                           --db-prefix {params.dbprefix} \
                           --input {input.rep} \
                           --output-prefix {output.bioboxes} \
                           --ranks {params.ranks} \
                           --output-format bioboxes {params.args} > {log} 2>&1
        mv {output.bioboxes}.tre {output.bioboxes}
        """