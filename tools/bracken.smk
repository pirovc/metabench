
rule bracken_build:
    output:
        db_kraken = temp("bracken/{vers}/{dtbs}/{dtbs_args}/database.kraken"),
        db_kraken_mers = temp("bracken/{vers}/{dtbs}/{dtbs_args}/database100mers.kraken"),
        kmer_dist = "bracken/{vers}/{dtbs}/{dtbs_args}/database100mers.kmer_distrib"
    benchmark:
        repeat("bracken/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    log:
        "bracken/{vers}/{dtbs}/{dtbs_args}.build.log"
    threads:
        config["threads"]
    conda:
        ("../envs/bracken.yaml")
    params:
        path = lambda wildcards: config["tools"]["bracken"][wildcards.vers],
        outprefix = "bracken/{vers}/{dtbs}/", #no args, to link from kraken2
        fixed_args = lambda wildcards: dict2args(config["run"]["bracken"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        kraken2db = lambda wildcards: config["run"]["bracken"][wildcards.vers][wildcards.dtbs]["fixed_args"]["-d"],
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        rm -rf "{params.outprefix}{wildcards.dtbs_args}" #  remove auto-generated folder
        ln -s "{params.kraken2db}" "{params.outprefix}" #  link kraken2 build
        {params.path}bracken-build -t {threads} {params.args} {params.fixed_args} > {log} 2>&1
        """

rule bracken_build_size:
    input:
        kmer_dist = "bracken/{vers}/{dtbs}/{dtbs_args}/database100mers.kmer_distrib"
    output:
        "bracken/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "{build_size_cmd} {input} > {output}"

rule bracken_profiling:
    input:
        rep = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep"
    output:
        bra = temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.bracken"),
        out = temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.out")
    log:
        "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    conda:
        ("../envs/bracken.yaml")
    params:
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/",
    shell:
        """
        bracken -d {params.dbprefix} -i {input.rep} -o {output.out} -w {output.bra} > {log} 2>&1
        """

rule bracken_profiling_format:
    input:
        bra = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.bracken"
    output:
        bioboxes = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes.log"
    conda:
        ("../envs/evals.yaml")
    params:
        scripts_path = os.path.join(workflow.basedir, "../scripts"),
        ranks = lambda wildcards: " ".join(config["default_ranks"]),
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/nodes.dmp",
                                            os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/names.dmp"],
        header = lambda wildcards: header_bioboxes_profiling("kraken2",
                                                           config["default_ranks"],
                                                           [os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/nodes.dmp",
                                                            os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/names.dmp"],
                                                           wildcards),
    shell: 
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}
        python3 {params.scripts_path}/profile.py \
                                    --taxid-col 4 \
                                    --perc-col 0 \
                                    --input-file {input.bra} \
                                    --taxonomy ncbi \
                                    --taxonomy-files {params.taxonomy_files} \
                                    --ranks >> {output.bioboxes} 2> {log}

        """
