import shutil
import os

snake_dir = workflow.basedir

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
        f"file:{snake_dir}/../wrappers/ganon/2.4.3/build"


rule ganon_build_size:
    input:
        multiext("ganon/{vers}/{dtbs}/{dtbs_args}/ganon_db", ".tax", ".hibf")
    output:
        "ganon/{vers}/{dtbs}/{dtbs_args}.build.size.tsv"
    shell:
        "{build_size_cmd} {input} > {output}"


# rule ganon_classify:
#     input:
#         paired_reads = lambda wildcards: [os.path.abspath(config["samples"][wildcards.samp]["fq1"]),
#                                           os.path.abspath(config["samples"][wildcards.samp]["fq2"])],
#         dbs = lambda wildcards: [os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.hibf')]
#     output:
#         rep = "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.rep",
#         one = temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.one"),
#         tre = temp("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.tre")
#     log:
#         "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.log"
#     benchmark:
#         repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}.binning.bench.tsv", config["repeat"])
#     threads:
#         config["threads"]
#     params:
#         extra = lambda wildcards: str2args(wildcards.b_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
#     wrapper:
#         f"file:{snake_dir}/../wrappers/ganon/latest/classify"

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
        batch_reads = temp(lambda wildcards: make_batch_reads_file(config["samples"], wildcards)),
        dbs = lambda wildcards: sorted([os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.hibf')])
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
    params:
        extra = lambda wildcards: str2args(wildcards.b_args) + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
    wrapper:
        f"file:{snake_dir}/../wrappers/ganon/2.4.3/classify"


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
        #db_tax = lambda wildcards: os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/ganon_db.tax"
        db_tax = lambda wildcards: [os.path.join(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/", f) for f in os.listdir(os.path.abspath(config["run"]["ganon"][wildcards.vers]["dbs"][wildcards.dtbs]) + "/" + wildcards.dtbs_args + "/") if f.endswith('.tax')]
    output:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bioboxes"
    log:
        "ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.log"
    benchmark:
        repeat("ganon/{vers}/{samp}/{dtbs}/{dtbs_args}/{b_args}/{p_args}.profiling.bench.tsv", config["repeat"])
    params:
        ranks = lambda wildcards: " ".join(config["default_ranks"]),
        extra = lambda wildcards: str2args(wildcards.p_args)# + " " + dict2args(config["run"]["ganon"][wildcards.vers]["fixed_args"])
    wrapper:
        f"file:{snake_dir}/../wrappers/ganon/2.4.3/report"
