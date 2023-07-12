rule kmcp_build:
    output:
        index_dir = directory("kmcp/{vers}/{dtbs}/{dtbs_args}/kmcp_db"),
        tmp_dir = temp(directory("kmcp/{vers}/{dtbs}/{dtbs_args}/kmcp_tmp")),
        taxid_map = "kmcp/{vers}/{dtbs}/{dtbs_args}/kmcp_db/taxid.map",
        name_map = "kmcp/{vers}/{dtbs}/{dtbs_args}/kmcp_db/name.map"
    benchmark:
        repeat("kmcp/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    log:
        "kmcp/{vers}/{dtbs}/{dtbs_args}.build.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/kmcp.yaml")
    params:
        path = lambda wildcards: config["tools"]["kmcp"][wildcards.vers],
        db = lambda wildcards: config["dbs"][wildcards.dtbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["kmcp"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        {params.path}kmcp compute \
                          --in-dir {params.db[folder]} \
                          --out-dir {output.tmp_dir} \
                          --threads {threads} \
                           {params.fixed_args} \
                           {params.args} > {log} 2>&1

        {params.path}kmcp index \
                          --in-dir {output.tmp_dir} \
                          --out-dir {output.index_dir} \
                          --threads {threads} >> {log} 2>&1
        
        # Make mapping files for KMCP targets "GCF_019263745.1_ASM1926374v1_genomic"
        find {params.db[folder]} -name "*.fna.gz" | xargs -I {{}} basename {{}} | awk '{{split($1,a,"_"); print a[1]"_"a[2]"\t"substr($1,1,length($1)-7)}}' > "{output.tmp_dir}/accver_filename.txt" 2>> {log}
        join <(sort -k 1,1 "{output.tmp_dir}/accver_filename.txt") <(sort -k 1,1 {params.db[assembly_summary]}) -t$'\t' -o "1.2,2.6" > {output.index_dir}/taxid.map 2>> {log}
        join <(sort -k 1,1 "{output.tmp_dir}/accver_filename.txt") <(sort -k 1,1 {params.db[assembly_summary]}) -t$'\t' -o "1.2,2.8" > {output.index_dir}/name.map 2>> {log}
        
        # Prepare taxonomy
        mkdir {output.index_dir}/taxonomy
        tar -xf {params.db[taxonomy_files]} -C {output.index_dir}/taxonomy/ >> {log} 2>&1
        """

rule kmcp_build_size:
    input:
        index_dir = "kmcp/{vers}/{dtbs}/{dtbs_args}/kmcp_db"
    output:
        "kmcp/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "du --bytes --dereference --max-depth 0 {input} > {output}"  # output in bytes


rule kmcp_binning:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output:
        search_out = temp("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tsv.gz")
    benchmark:
        repeat("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/kmcp.yaml")
    params:
        path = lambda wildcards: config["tools"]["kmcp"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["kmcp"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args)
    shell:
        """
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}kmcp search \
                              --db-dir {params.dbprefix} \
                              --out-file {output.search_out} \
                              --threads {threads} \
                              {params.fixed_args} \
                              {params.args} \
                              {input.fq1} > {log} 2>&1
        else # paired-end
            {params.path}kmcp search \
                              --db-dir {params.dbprefix} \
                              --out-file {output.search_out} \
                              --threads {threads} \
                              {params.fixed_args} \
                              {params.args} \
                              --read1 {input.fq1} --read2 {params.input_fq2} > {log} 2>&1
        fi

        """

rule kmcp_binning_format:
    input: 
        search_out = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tsv.gz",
        taxid_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxid.map"
    output:
        bioboxes = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
    params:
        header = lambda wildcards: header_bioboxes_binning("kmcp", wildcards)
    shell:
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # Get taxid from target from taxid.map file
        # search_out fields: #query qLen    qKmers  FPR hits    target  chunkIdx    chunks  tLen    kSize   mKmers  qCov    tCov    jacc    queryIdx
        # Check if end of read id is "/1" and remove it

        join -1 2 -2 1 -t$'\t' -o "1.1,2.2,1.2" \
             <(zcat {input.search_out} | tail -n+2 | cut -f 1,6 | sort -k 2,2) \
             <(sort -k 1,1 {input.taxid_map})  | \
             awk '{{split($3,a,"_");
                    accver=a[1]"_"a[2];
                    if(substr($1,length($1)-1)=="/1"){{
                        $1=substr($1,0,length($1)-2);
                    }};
                    print $1"\t"$2"\t"accver}}' >> {output.bioboxes}
        """


rule kmcp_profiling:
    input:
        search_out="kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tsv.gz",
        taxid_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxid.map",
        name_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/name.map"
    output:
        bioboxes = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    benchmark:
        repeat("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    log:
        "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    threads:
        config["threads"] if config["threads"] <= 4 else 4
    conda:
        srcdir("../envs/kmcp.yaml")
    params:
        path = lambda wildcards: config["tools"]["kmcp"][wildcards.vers],
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxonomy"],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db",
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        kmcp profile {input.search_out} \
                     --threads {threads} \
                     --cami-report {output.bioboxes} \
                     --taxid-map {input.taxid_map}  \
                     --name-map {input.name_map} \
                     --taxdump {params.taxonomy_files} {params.args} > {log} 2>&1
        # remove .profile suffix
        mv {output.bioboxes}.profile {output.bioboxes}
        """