rule kraken2_build:
    output: db1="kraken2/{vers}/{dtbs}/{prms}/hash.k2d",
            db2="kraken2/{vers}/{dtbs}/{prms}/opts.k2d",
            db3="kraken2/{vers}/{dtbs}/{prms}/taxo.k2d",
            db4="kraken2/{vers}/{dtbs}/{prms}/taxonomy/nodes.dmp",
            db5="kraken2/{vers}/{dtbs}/{prms}/taxonomy/names.dmp",
            fasta=temp("kraken2/{vers}/{dtbs}/{prms}/input.fasta"),
            accession2taxid=temp("kraken2/{vers}/{dtbs}/{prms}/taxonomy/kraken2.accession2taxid")
    benchmark: "kraken2/{vers}/{dtbs}/{prms}.build.bench.tsv"
    log: "kraken2/{vers}/{dtbs}/{prms}.log"
    threads: config["threads"]
    conda: srcdir("../envs/kraken2.yaml") 
    params: outprefix = "kraken2/{vers}/{dtbs}/{prms}/",
            path = lambda wildcards: config["tools"]["kraken2"][wildcards.vers],
            folder = lambda wildcards: os.path.abspath(config["dbs"][wildcards.dtbs]["folder"]),
            extens = lambda wildcards: config["dbs"][wildcards.dtbs]["extension"],
            detail = lambda wildcards: config["dbs"][wildcards.dtbs]["details"],
            taxonomy_files = " ".join(config["taxonomy_files"]) if "taxonomy_files" in config else "",
            fixed_params = lambda wildcards: config["run"]["kraken2"]["fixed_params"] if "fixed_params" in config["run"]["kraken2"] else "",
            params = lambda wildcards: str2params(wildcards.prms),
    shell: 
        """
        mkdir -p "{params.outprefix}taxonomy"
        mkdir -p "{params.outprefix}library"
        tar xf {params.taxonomy_files} -C {params.outprefix}taxonomy/ nodes.dmp names.dmp > {log} 2>&1
        find {params.folder} -name *{params.extens} | xargs zcat > {output.fasta} 2>> {log}
        awk 'BEGIN {{FS="\\t";OFS="\\t"; print "accession","accession.version","taxid","gi"}}{{ split($4,acc,"."); print acc[1],$4,$6,"0" }}' {params.detail} > {output.accession2taxid} 2>> {log}
        {params.path}kraken2-build --db {params.outprefix} --no-masking --add-to-library {output.fasta} >> {log} 2>&1
        {params.path}kraken2-build --build --db {params.outprefix} --threads {threads} {params.params} {params.fixed_params} >> {log} 2>&1
        rm -rfv {params.outprefix}library/ {params.outprefix}seqid2taxid.map {params.outprefix}taxonomy/prelim_map.txt >> {log} 2>&1
        """

rule kraken2_build_size:
    input: "kraken2/{vers}/{dtbs}/{prms}/hash.k2d",
           "kraken2/{vers}/{dtbs}/{prms}/opts.k2d",
           "kraken2/{vers}/{dtbs}/{prms}/taxo.k2d",
           "kraken2/{vers}/{dtbs}/{prms}/taxonomy/nodes.dmp",
           "kraken2/{vers}/{dtbs}/{prms}/taxonomy/names.dmp"
    output: "kraken2/{vers}/{dtbs}/{prms}.size.tsv"
    shell: "du --block-size=1 {input} > {output}"  # output in bytes

rule kraken2_classify:
    input: fq1 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
           dbfolder = lambda wildcards: os.path.abspath(config["run"]["kraken2"]["dbs"][wildcards.dtbs])
    output: res=temp("kraken2/{vers}/{samp}/{dtbs}/{prms}.res")
    benchmark: "kraken2/{vers}/{samp}/{dtbs}/{prms}.classify.bench.tsv"
    log: "kraken2/{vers}/{samp}/{dtbs}/{prms}.log"
    threads: config["threads"]
    conda: srcdir("../envs/kraken2.yaml")
    params: outprefix = "kraken2/{vers}/{samp}/{dtbs}/{prms}",
            path = lambda wildcards: config["tools"]["kraken2"][wildcards.vers],
            input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
            fixed_params = lambda wildcards: config["run"]["kraken2"]["fixed_params"] if "fixed_params" in config["run"]["kraken2"] else "",
            params = lambda wildcards: str2params(wildcards.prms)
    shell:
        """
        if [[ ! -z "{params.path}" ]]; then
            source deactivate;
        fi
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            {params.path}kraken2 --db {input.dbfolder} --threads {threads} --output {output.res} --gzip-compressed {params.params} {input.fq1} > {log} 2>&1
        else # paired-end
            {params.path}kraken2 --db {input.dbfolder} --threads {threads} --output {output.res} --gzip-compressed {params.params} --paired {input.fq1} {params.input_fq2} > {log} 2>&1
        fi
        """

rule kraken2_classify_format:
    input: res="kraken2/{vers}/{samp}/{dtbs}/{prms}.res",
    output: "kraken2/{vers}/{samp}/{dtbs}/{prms}.bioboxes"
    params: input_fq2 = lambda wildcards: os.path.abspath(config["samples"][wildcards.samp]["fq2"]) if config["samples"][wildcards.samp]["fq2"] else "",
    shell:
        """
        # bioboxes header
        printf "@Version:0.9.1\\n@SampleID:kraken2 {wildcards.vers} {wildcards.samp} {wildcards.dtbs} {wildcards.prms}\\n@@SEQUENCEID\\tBINID\\tTAXID\\n" > {output}

        # output header changes when sinle or paired (/1 or nothing)
        if [[ -z "{params.input_fq2}" ]]; then # single-end
            header_suffix=2;
        else # paired-end
            header_suffix=0;
        fi

        grep "^C" {input.res} | awk -v header_suffix="${{header_suffix}}" 'FS="\\t"{{print substr($2,1,length($2)-header_suffix)"\\t0\\t"$3}}' >> {output}
        """
