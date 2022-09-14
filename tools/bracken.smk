rule bracken_build:
    output:
        "bracken/{vers}/{dtbs}/{args}/database.kraken"
    benchmark:
        repeat("bracken/{vers}/{dtbs}/{args}.build.bench.tsv", config["repeat"])
    log:
        "bracken/{vers}/{dtbs}/{args}.build.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/bracken.yaml")
    params:
        path = lambda wildcards: config["tools"]["bracken"][wildcards.vers],
        outprefix = "bracken/{vers}/{dtbs}/{args}/",
        fixed_args = lambda wildcards: dict2args(config["run"]["bracken"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        kraken2db = lambda wildcards: config["run"]["bracken"][wildcards.vers][wildcards.dtbs]["fixed_args"]["-d"],
        args = lambda wildcards: str2args(wildcards.args)
    shell: 
        """
        ln -s {params.kraken2db}/* {params.outprefix}
        {params.path}bracken-build -t {threads} {params.args} {params.fixed_args} > {log} 2>&1
        """

rule bracken_build_size:
    input:
        "bracken/{vers}/{dtbs}/{args}/database.kraken"
    output:
        "bracken/{vers}/{dtbs}/{args}.build.size.tsv"
    shell:
        "du --block-size=1 --dereference {input0} > {output}"  # output in bytes


rule bracken_profile:
    input:
        rep = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.rep"
    output:
        bra = temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.bracken"),
        out = temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.out")
    log:
        "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.log"
    conda:
        srcdir("../envs/bracken.yaml")
    params:
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/",
    shell:
        """
        bracken -d {params.dbprefix} -i {input.rep} -o {output.out} -w {output.bra} > {log} 2>&1
        """

rule bracken_profile_format:
    input:
        bra = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.bracken"
    output:
        bioboxes = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.bioboxes"
    log:
        "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{args}.profile.bioboxes.log"
    conda:
        srcdir("../envs/evals.yaml")
    params:
        scripts_path = srcdir("../scripts/"),
        ranks = " ".join(config["profile_ranks"]),
        taxonomy_files = lambda wildcards: [os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/nodes.dmp",
                                            os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/names.dmp"],
        header = lambda wildcards: header_bioboxes_profile("kraken2",
                                                           config["profile_ranks"],
                                                           [os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/nodes.dmp",
                                                            os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/taxonomy/names.dmp"],
                                                           wildcards),
    shell: 
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}
        python3 {params.scripts_path}profile.py \
                                    --taxid-col 4 \
                                    --perc-col 0 \
                                    --input-file {input.bra} \
                                    --taxonomy ncbi \
                                    --taxonomy-files {params.taxonomy_files} \
                                    --ranks {params.ranks} >> {output.bioboxes} 2> {log}

        """