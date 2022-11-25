rule kraken2_build:
    output:
        db1="kraken2/{vers}/{dtbs}/{dtbs_args}/hash.k2d",
        db2="kraken2/{vers}/{dtbs}/{dtbs_args}/opts.k2d",
        db3="kraken2/{vers}/{dtbs}/{dtbs_args}/taxo.k2d",
        db4="kraken2/{vers}/{dtbs}/{dtbs_args}/taxonomy/nodes.dmp",
        db5="kraken2/{vers}/{dtbs}/{dtbs_args}/taxonomy/names.dmp",
        fasta=temp("kraken2/{vers}/{dtbs}/{dtbs_args}/input.fasta"),
        accession2taxid=temp("kraken2/{vers}/{dtbs}/{dtbs_args}/taxonomy/kraken2.accession2taxid")
    benchmark:
        repeat("kraken2/{vers}/{dtbs}/{dtbs_args}.build.bench.tsv", config["repeat"])
    log:
        "kraken2/{vers}/{dtbs}/{dtbs_args}.build.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/kraken2.yaml")
    priority: 1  # to run before bracken
    params:
        path = lambda wildcards: config["tools"]["kraken2"][wildcards.vers],
        outprefix = "kraken2/{vers}/{dtbs}/{dtbs_args}/",
        db = lambda wildcards: config["dbs"][wildcards.dtbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["kraken2"][wildcards.vers][wildcards.dtbs]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.dtbs_args)
    shell: 
        """
        #mkdir -p "{params.outprefix}taxonomy"
        #mkdir -p "{params.outprefix}library"
        #tar xf {params.db[taxonomy_files]} -C {params.outprefix}taxonomy/ nodes.dmp names.dmp > {log} 2>&1
        #awk 'BEGIN {{FS="\\t";OFS="\\t"; print "accession","accession.version","taxid","gi"}}{{ split($4,acc,"."); print acc[1],$4,$6,"0" }}' {params.db[details]} > {output.accession2taxid} 2>> {log}
        # Kraken2
        find {params.db[folder]} -name *{params.db[extension]} | xargs zcat > {output.fasta} 2>> {log}
        {params.path}kraken2-build --db {params.outprefix} --no-masking --add-to-library {output.fasta} >> {log} 2>&1
        {params.path}kraken2-build --db {params.outprefix} --download-taxonomy >> {log} 2>&1
        {params.path}kraken2-build --build --db {params.outprefix} --threads {threads} {params.args} {params.fixed_args} >> {log} 2>&1
        rm -rfv {params.outprefix}taxonomy/prelim_map.txt >> {log} 2>&1
        """

rule kraken2_build_size:
    input:
        "kraken2/{vers}/{dtbs}/{dtbs_args}/hash.k2d",
        "kraken2/{vers}/{dtbs}/{dtbs_args}/opts.k2d",
        "kraken2/{vers}/{dtbs}/{dtbs_args}/taxo.k2d",
        "kraken2/{vers}/{dtbs}/{dtbs_args}/taxonomy/nodes.dmp",
        "kraken2/{vers}/{dtbs}/{dtbs_args}/taxonomy/names.dmp"
    output:
        "kraken2/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "du --block-size=1 {input} > {output}"  # output in bytes

rule kraken2_classify:
    input:
        fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"])
    output:
        res=temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.res"),
        rep=temp("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep")
    benchmark:
        repeat("kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
    log:
        "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
    threads:
        config["threads"]
    conda:
        srcdir("../envs/kraken2.yaml")
    params:
        path = lambda wildcards: config["tools"]["kraken2"][wildcards.vers],
        dbprefix = lambda wildcards: os.path.abspath(config["run"]["kraken2"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/",
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        fixed_args = lambda wildcards: dict2args(config["run"]["kraken2"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args)
    shell:
        """
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}kraken2 --db {params.dbprefix} \
                                 --output {output.res} \
                                 --report {output.rep} \
                                 --threads {threads} \
                                 --gzip-compressed \
                                 {params.args} {params.fixed_args} \
                                 {input.fq1} > {log} 2>&1
        else # paired-end
            {params.path}kraken2 --db {params.dbprefix} \
                                 --output {output.res} \
                                 --report {output.rep} \
                                 --threads {threads} \
                                 --gzip-compressed \
                                 {params.args} {params.fixed_args} \
                                 --paired {input.fq1} {params.input_fq2} > {log} 2>&1
        fi
        """

rule kraken2_classify_format:
    input:
        res = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.res",
    output:
        bioboxes = "kraken2/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bioboxes"
    params:
        input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
        header = lambda wildcards: header_bioboxes_binning("kraken2", wildcards)
    shell:
        """
        # bioboxes header
        echo "{params.header}" > {output.bioboxes}

        # output header changes when sinle or paired (/1 or nothing)
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            header_suffix=2;
        else # paired-end
            header_suffix=0;
        fi

        grep "^C" {input.res} | awk -v header_suffix="${{header_suffix}}" 'FS="\\t"{{print substr($2,1,length($2)-header_suffix)"\\t"$3}}' >> {output.bioboxes}
        """

# running with bracken
# rule kraken2_profiling_format:
