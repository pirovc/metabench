rule metacache_build:
    output:
        db1 = "metacache/{vers}/{dtbs}/{dtbs_args}/metacache_db.cache0",
        db2 = "metacache/{vers}/{dtbs}/{dtbs_args}/metacache_db.meta",
        tmp_dir = temp(directory("metacache/{vers}/{dtbs}/{dtbs_args}/metacache_tmp/")),
    benchmark:
        repeat("metacache/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    log:
        "metacache/{vers}/{dtbs}/{dtbs_args}.build.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/metacache.yaml")
    params:
        path = lambda wildcards: config["tools"]["metacache"][wildcards.vers],
        outprefix = "metacache/{vers}/{dtbs}/{dtbs_args}/metacache_db",
        db = lambda wildcards: config["dbs"][wildcards.dtbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["metacache"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        # Tmp folder for taxonomy files
        mkdir -p {output.tmp_dir}
        tar -xf {params.db[taxonomy_files]} -C {output.tmp_dir} > {log} 2>&1

        # Download acc2txid
        {params.path}download-ncbi-taxmaps {output.tmp_dir} >> {log} 2>&1

        {params.path}metacache build {params.outprefix} \
                                     {params.db[folder]} \
                                     -taxonomy {output.tmp_dir} \
                                     -taxpostmap {output.tmp_dir}/*.accession2taxid \
                                     {params.fixed_args} \
                                     {params.args} >> {log} 2>&1
        """

rule metacache_build_size:
    input:
        db1 = "metacache/{vers}/{dtbs}/{dtbs_args}/metacache_db.cache0",
        db2 = "metacache/{vers}/{dtbs}/{dtbs_args}/metacache_db.meta"
    output:
        "metacache/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "du --bytes --dereference {input} > {output}"  # output in bytes


rule metacache_classify:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output:
        query_out = temp("metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.out")
    benchmark:
        repeat("metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/metacache.yaml")
    params:
        path = lambda wildcards: config["tools"]["metacache"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["metacache"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/metacache_db",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["metacache"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args)
    shell:
        """
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}metacache query {params.dbprefix} \
                                         {input.fq1} \
                                         -no-summary -no-query-params \
                                         -taxids -mapped-only -separate-cols -separator '\t' \
                                         -threads {threads} {params.fixed_args} {params.args} > {output.query_out} 2> {log}
        else # paired-end
            {params.path}metacache query {params.dbprefix} \
                                         {input.fq1} {params.input_fq2} -pairfiles \
                                         -no-summary -no-query-params \
                                         -taxids -mapped-only -separate-cols -separator '\t' \
                                         -threads {threads} {params.fixed_args} {params.args} > {output.query_out} 2> {log}
        fi
        """

rule metacache_classify_format:
    input: 
        query_out = "metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.out"
    output:
        bioboxes = "metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes",
        taxid_map = temp("metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.map")
    conda:
        srcdir("../envs/metacache.yaml")
    params:
        path = lambda wildcards: config["tools"]["metacache"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["metacache"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/metacache_db",
        header = lambda wildcards: header_bioboxes_binning("metacache", wildcards)
    shell:
        """

        # Get DB info to link taxids and sequence ids
        # cols: name    sequence    form    variety subspecies  species subgenus    genus   subtribe    tribe   subfamily   family  suborder    ordersubclass   class   subphylum   phylum  subkingdom  kingdom domain
        {params.path}metacache info {params.dbprefix} lin 2> /dev/null | tail -n+2 | cut -f 2,6 > {output.taxid_map}

        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # Print taxid classification (not sequence) - check if /1 on readid and remove it
        tail -n+3 {input.query_out} | awk 'BEGIN{{FS=OFS="\t"}}{{if($2!="sequence"){{if(substr($1,length($1)-1)=="/1"){{$1=substr($1,0,length($1)-2);}}; print $1,$4}}}}' >> {output.bioboxes}

        # Print sequence classification (get taxid from map) - check if /1 on readid and remove it
        # out: readid <tab> taxid <tab> <tab> sequence id (2.99 placeholder to insert an empty field)
        join -1 4 -2 1 <(tail -n+3 {input.query_out} | awk 'BEGIN{{FS=OFS="\t"}}{{if($2=="sequence"){{print}}}}' | sort -k 4,4 -t$'\t') <(sort -k 1,1 -t$'\t' {output.taxid_map} ) -o "1.1,2.2,2.99,1.3" -t$'\t' >> {output.bioboxes}
        """


rule metacache_profiling:
    input:
        bioboxes = "metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
    output:
        bioboxes = touch("metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes")
    benchmark:
        repeat("metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    log:
        "metacache/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    shell:
        """
        """