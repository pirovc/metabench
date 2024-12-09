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
        ("../envs/kmcp.yaml")
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
                          --force \
                           {params.fixed_args} \
                           {params.args} > {log} 2>&1

        {params.path}kmcp index \
                          --in-dir {output.tmp_dir} \
                          --out-dir {output.index_dir} \
                          --force \
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
        "{build_size_cmd} {input} > {output}"


rule kmcp_binning:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
        taxid_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxid.map",
        name_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/name.map"
    output:
        search_out = temp("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tsv.gz"),
        b_bioboxes = temp("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.gz")
    benchmark:
        repeat("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    threads:
        config["threads"]
    conda:
        ("../envs/kmcp.yaml")
    params:
        path = lambda wildcards: config["tools"]["kmcp"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db",
        outprefix = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxonomy"],
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
    
        # kmcp adds .binning.gz to --binning-result
        kmcp profile {output.search_out} \
                     --threads {threads} \
                     --binning-result {params.outprefix} \
                     --taxid-map {input.taxid_map}  \
                     --name-map {input.name_map} \
                     --taxdump {params.taxonomy_files} >> {log} 2>&1
        """

rule kmcp_profiling:
    input:
        search_out = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tsv.gz",
        taxid_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxid.map",
        name_map = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/name.map"
    output:
        p_bioboxes = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes",
    benchmark:
        repeat("kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    log:
        "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    threads:
        # as stated in the log of kmcp profile:
        # 10:10:55.163 [INFO] using a lot of threads does not always accelerate processing, 4-threads is fast enough
        # 10:10:55.163 [INFO] kmcp v0.9.2
        # 10:10:55.163 [INFO]   https://github.com/shenwei356/kmcp
        config["threads"] if config["threads"] <= 4 else 4
    conda:
        ("../envs/kmcp.yaml")
    params:
        path = lambda wildcards: config["tools"]["kmcp"][wildcards.vers],
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db/taxonomy"],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kmcp"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/kmcp_db",
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        kmcp profile {input.search_out} \
                     --threads {threads} \
                     --cami-report {output.p_bioboxes} \
                     --taxid-map {input.taxid_map}  \
                     --name-map {input.name_map} \
                     --taxdump {params.taxonomy_files} {params.args} > {log} 2>&1
        
        # remove .profile suffix
        mv {output.p_bioboxes}.profile {output.p_bioboxes}
        """

rule kmcp_binning_format:
    input:
        bioboxes = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.gz"
    output:
        bioboxes = "kmcp/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
    params:
        header = lambda wildcards: header_bioboxes_binning("kmcp", wildcards)
    shell:
        """
        # "catch exit status 1" grep wrapper
        c1zgrep() {{ zgrep "$@" || test $? = 1; }}

        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # Check if end of read id is "/1" and remove it
        c1zgrep -e "^@" -e "^#" -v {input.bioboxes} | awk 'BEGIN{{FS=OFS="\t"}}
            {{
            if(substr($1,length($1)-1)=="/1"){{
                $1=substr($1,0,length($1)-2);
            }};

            print $1,$2,"","";
            }}' >> {output.bioboxes}
        """
