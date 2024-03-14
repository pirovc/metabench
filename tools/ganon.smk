rule ganon_build:
    input:
        lambda wildcards: config["dbs"][wildcards.dtbs]["folder"]
    output:
        multiext("ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db", ".tax", ".hibf")
    log:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.log"
    benchmark:
        repeat("ganon/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    threads:
        config["threads"]
    params:
        extension = lambda wildcards: config["dbs"][wildcards.dtbs]["extension"],
        taxonomy = lambda wildcards: config["dbs"][wildcards.dtbs]["taxonomy"],
        taxonomy_files = lambda wildcards: config["dbs"][wildcards.dtbs]["taxonomy_files"],
        ncbi_file_info = lambda wildcards: config["dbs"][wildcards.dtbs]["assembly_summary"],
        extra = lambda wildcards: str2args(wildcards.dtbs_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers][wildcards.dtbs]["fixed_args"])
    wrapper:
        "file:///home/pirovc/code/metabench/wrappers/ganon/2.1.0/build"


rule ganon_build_size:
    input:
        multiext("ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db", ".tax", ".hibf")
    output:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "{build_size_cmd} {input} > {output}"


rule ganon_classify:
    input:
        paired_reads = lambda wildcards: [os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
                                          os.path.abspath(config["samples"][wildcards.samp]["fq2"])],
        db = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.hibf"
    output:
        rep = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
        one = temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one"),
        tre = temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tre")
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    threads:
        config["threads"]
    params:
        extra = lambda wildcards: str2args(wildcards.b_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
    wrapper:
        "file:///home/pirovc/code/metabench/wrappers/ganon/2.1.0/classify"

rule ganon_binning_format:
    input: 
        one = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one",
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

rule ganon_report:
    input:
        rep = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
        db_tax = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"
    output:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    params:
        ranks = lambda wildcards: " ".join(config["default_ranks"]),
        extra = lambda wildcards: str2args(wildcards.p_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
    wrapper:
        "file:///home/pirovc/code/metabench/wrappers/ganon/2.1.0/report"
