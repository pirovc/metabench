import shutil
import os

rule git:
    output: "{commit}/ganon/README.md"
    log: git_log = os.path.abspath("{commit}/git.log")
    shell:
        """
        cd {wildcards.commit}
        git clone https://github.com/pirovc/ganon.git > {log} 2>&1
        cd ganon
        git reset --hard {wildcards.commit} >> {log} 2>&1
        git submodule update --init --recursive >> {log} 2>&1

        conda create --yes -c bioconda -c conda-forge -n "env_ganon_{wildcards.commit}" raptor pip "gcc==14.4.0" cxx-compiler zlib bzip2 "cmake>=3.5" make >> {log} 2>&1
        """

rule comp:
    input: "{commit}/ganon/README.md"
    output: ganon_classify_exec = "{commit}/ganon/build/ganon-classify",
            ganon_build_exec = "{commit}/ganon/build/ganon-build"
    log: os.path.abspath("{commit}/comp.log")
    conda:
        ("../envs/ganon.yaml")
    params: build_folder = "{commit}/ganon/build",
            #compile_params = lambda wildcards: config["commits"][wildcards.commit]["compile_params"]
    shell: 
        """
        conda activate "env_ganon_{wildcards.commit}"
        mkdir -p {params.build_folder}
        cd {params.build_folder}
        cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. > {log} 2>&1
        make -j >> {log} 2>&1
        make install >> {log} 2>&1
        cd ..
        pip install --pre .[dev]  >> {log} 2>&1
        """

rule test:
    input: ganon_classify_exec = "{commit}/ganon/build/ganon-classify",
            ganon_build_exec = "{commit}/ganon/build/ganon-build"
    output: env = "env_ganon_{commit}"
    log: os.path.abspath("{commit}/test.log")
    conda:
        ("../envs/ganon.yaml")
    params: build_folder = "{commit}/ganon/build"
    shell: 
        """
        base_dir="$(pwd)"
        conda activate "env_ganon_{wildcards.commit}"
        cd {params.build_folder}
        ctest -VV . > {log} 2>&1
        cd ..
        #python3 -m unittest discover -s tests/ganon/integration/ >> {log} 2>&1

        cd "${{base_dir}}"
        ganon -v > {output.env} 2>&1
        ganon-build -v >> {output.env} 2>&1
        ganon-classify -v >> {output.env} 2>&1
        conda list >> {output.env} 2>&1

        # Remove code
        #rm -r {wildcards.commit}
        """

def define_env(ganon_cfg):
    return "env_ganon_" + ganon_cfg if ganon_cfg else []

def make_batch_reads_file(samples, wildcards):
    batch_reads_file = f"{wildcards.vers}_{wildcards.dtbs}_{wildcards.dtbs_args}_{wildcards.b_args}.tsv"
    with open(batch_reads_file, "w") as brf:
        for sample, f in samples.items():
            fq1s = f["fq1"].split(" ")
            fq2s = f["fq2"].split(" ")
            for i, fq1 in enumerate(fq1s):
                fq2 = fq2s[i] if fq2s != [""] else ""
                print(f"ganon/{wildcards.vers}/{sample}/{wildcards.dtbs}/{wildcards.dtbs_args}/{wildcards.b_args}", fq1, fq2, sep="\t", file=brf)
    return batch_reads_file

rule ganon_classify_batch:
    input:
        batch_reads = lambda wildcards: make_batch_reads_file(config["samples"], wildcards),
        dbs = lambda wildcards: sorted([os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.hibf')]),
        env = lambda wildcards: define_env(config["tools"]["ganon"][wildcards.vers])
    output:
        rep = ["ganon/{vers}/" + samp + "/{dtbs}/{dtbs_args}/{b_args}.rep" for samp in config["samples"]],
        one = temp(["ganon/{vers}/" + samp + "/{dtbs}/{dtbs_args}/{b_args}.one" for samp in config["samples"]]),
        tre = temp(["ganon/{vers}/" + samp + "/{dtbs}/{dtbs_args}/{b_args}.tre" for samp in config["samples"]])
    log:
        "ganon_{vers}_{dtbs}_{dtbs_args}_{b_args}.binning.log"
    benchmark:
       repeat("ganon_{vers}_{dtbs}_{dtbs_args}_{b_args}.binning.bench.tsv", config["repeat"])
    threads:
        config["threads"]
    conda:
        ("../envs/ganon.yaml")
    params:
        dbprefixes = lambda wildcards, input: [db[:-5] for db in input.dbs],
        fixed_args = lambda wildcards: dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"]),
        args = lambda wildcards: str2args(wildcards.b_args)
    shell:
        """
        if [[ -n {input.env} ]]; then
            conda activate {input.env}
        fi

        ganon classify \
                --db-prefix {params.dbprefixes} \
                --batch-reads {input.batch_reads} \
                --output-prefix ./ \
                --threads {threads} \
                --output-one \
                --verbose \
                {params.fixed_args} \
                {params.args} > {log} 2>&1
        
        rm {input.batch_reads}
        """

rule mv_log_bench_batch:
    input:
        batch_log = "ganon_{vers}_{dtbs}_{dtbs_args}_{b_args}.binning.log",
        batch_bench = "ganon_{vers}_{dtbs}_{dtbs_args}_{b_args}.binning.bench.tsv"
    output:
        benchs = (["ganon/{vers}/" + samp + "/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv" for samp in config["samples"]]),
        logs = (["ganon/{vers}/" + samp + "/{dtbs}/{dtbs_args}/{b_args}.binning.log" for samp in config["samples"]])
    run:
        for bench in output.benchs:
            shutil.copyfile(input.batch_bench, bench)
        for lg in output.logs:
            shutil.copyfile(input.batch_log, lg)
        os.remove(input.batch_bench)
        os.remove(input.batch_log)
        
rule ganon_binning_format:
    input: 
        one = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one",
        dbtax = lambda wildcards: [os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.tax')]
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
        <(cut -f 1,2,3 {input.dbtax} | sort | uniq | sort -t$'\t' -k 1,1 {input.dbtax}) | \
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
        db_tax = lambda wildcards: [os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.tax')],
        env = lambda wildcards: define_env(config["tools"]["ganon"][wildcards.vers])
    output:
        bioboxes = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    conda:
        ("../envs/ganon.yaml")
    params:
        ranks = lambda wildcards: " ".join(config["default_ranks"]),
        args = lambda wildcards: str2args(wildcards.p_args)
    shell:
        """
        if [[ -n {input.env} ]]; then
            conda activate {input.env}
        fi
        ganon report \
            --db-prefix {input.db_tax} \
            --input {input.rep} \
            --output-prefix {output.bioboxes} \
            --ranks {params.ranks} \
            --output-format bioboxes {params.args} > {log} 2>&1
        mv {output.bioboxes}.tre {output.bioboxes}
        """