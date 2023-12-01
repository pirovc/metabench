# MetaBench

MetaBench is a pipeline to continuously benchmark metagenomics analysis tools. It covers database construction (*build*), taxonomic *binning* and *profiling*.

It supports:
- multiple tools
	- in several release versions
- multiple databases
	- with optional range of parameters for building
- multiple samples
	- for binning and profiling
	- both with optional range of parameters
- multiple evaluation metrics and threshold
- performance benchmarks with optional repeats

It outputs:
- standardized JSON files for integration
- [bioboxes](https://github.com/bioboxes/rfc/tree/master/data-format) output files for binning and profiling
- interactive dashboard to analyze results

It requires:
- Reference sequences to build databases
- Set of reads and ground truth in the bioboxes format to evaluate

Current configured tools:
- ganon
- kmcp
- metacache
- kraken2
- bracken

## Installation and requirements

MetaBench is written in Snakemake and makes use of conda/mamba internally to install dependencies. It uses Bokeh to plot the interactive dashboard.

```sh
pip install randomname
mamba install snakemake pandas
mamba install "bokeh==2.4.3"
git clone https://github.com/pirovc/metabench.git
cd metabench
```
## Usage example

### Build

Downloading some references for the build with [genome_updater](https://github.com/pirovc/genome_updater):

```sh
genome_updater.sh -d refseq -g bacteria -c "reference genome" -f "genomic.fna.gz" -o bac_rs -b refgen -t 8 -a
```

Create `config/build_test.yaml`:

```yaml
workdir: "build_example/"
threads: 8
repeat: 1

tools:
  ganon:
    "latest": ""
  kmcp:
    "latest": ""

dbs:
  "bac_rs_refgen":
    folder: "../bac_rs/refgen/files/"
    extension: ".fna.gz"
    taxonomy: "ncbi"
    taxonomy_files: "../bac_rs/refgen/taxdump.tar.gz"
    assembly_summary: "../bac_rs/refgen/assembly_summary.txt"

run:
   ganon:
     "latest":
       bac_rs_refgen:
         fixed_args:
           "--ncbi-file-info": "../bac_rs/refgen/assembly_summary.txt"
         args:
           "--max-fp": [0.0001, ""]
   kmcp:
     "latest":
       bac_rs_refgen:
         fixed_args:
         args:
```

In the example above MetaBench is set to build databases for 2 tools (ganon and kmcp). kmcp will run with default parameters only (no `args:`) and ganon will run with `--max-fp 0.0001` and default parameters.

Verify run with `--dry-run`:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda`

If everything finished correctly, the following files will be created (`tree -A build_test`):

```
build_test/
├── ganon
│   └── latest
│       └── bac_rs_refgen
│           ├── default
│           │   ├── ganon_db.hibf
│           │   ├── ganon_db.ibf -> ganon_db.hibf
│           │   └── ganon_db.tax
│           ├── default.build.bench.json
│           ├── default.build.bench.tsv
│           ├── default.build.log
│           ├── default.build.size.tsv
│           ├── --max-fp=0.0001
│           │   ├── ganon_db.hibf
│           │   ├── ganon_db.ibf -> ganon_db.hibf
│           │   └── ganon_db.tax
│           ├── --max-fp=0.0001.build.bench.json
│           ├── --max-fp=0.0001.build.bench.tsv
│           ├── --max-fp=0.0001.build.log
│           └── --max-fp=0.0001.build.size.tsv
└── kmcp
    └── latest
        └── bac_rs_refgen
            ├── default
            │   └── kmcp_db
            │       ├── name.map
            │       ├── R001
            │       │   ├── _block001.uniki
            │       │   ├── _block002.uniki
            │       │   ├── __db.yml
            │       │   └── __name_mapping.tsv
            │       ├── taxid.map
            │       └── taxonomy
            │           ├── citations.dmp
            │           ├── delnodes.dmp
            │           ├── division.dmp
            │           ├── gc.prt
            │           ├── gencode.dmp
            │           ├── images.dmp
            │           ├── merged.dmp
            │           ├── names.dmp
            │           ├── nodes.dmp
            │           └── readme.txt
            ├── default.build.bench.json
            ├── default.build.bench.tsv
            ├── default.build.log
            └── default.build.size.tsv
```

- `*.build.bench.json` contains the standardized metrics in JSON format. If `repeat > 1` in the config file, only the fastest run is selected.
- `*.build.bench.tsv` contains the raw benchmark metrics from Snakemake. If `repeat > 1` in the config file, one line for each run will be reported.
- `*.build.log` contains the STDOUT and STDERR from the run.
- `*.build.size.tsv` contains the size in bytes for the mandatory database files (`du --bytes`).

Obs: note that if no arguments are used in `args:` section of the configuration, the database folder/files will be named `default`. If parameters are used, databases are created based on them (`--max-fp 0.0001` -> `--max-fp=0.0001`, if more than one, connected by underscore "\_"). Any information provided in `fixed_args:` is not accounted for file/folder names.

Check the `config/build_example.yaml` for more examples on how to use the configuration file. Multiple databases, range of parameters and others can be configured to be executed in the same run.

### Classify (binning + profiling)

Classification includes both binning and profiling procedures. It requires databases (as created in the build process above) and one or more samples with single or paired `fastq` files.


Downloading some references for the build with [genome_updater](https://github.com/pirovc/genome_updater):

```sh
genome_updater.sh -d refseq -g bacteria -c "reference genome" -f "genomic.fna.gz" -o bac_rs -b refgen -t 8 -a
```

Create `config/classify_test.yaml`:

```yaml
workdir: "classify_test/"
threads: 8
repeat: 1

tools:
  ganon:
    "latest": ""
  kmcp:
    "latest": ""

samples:
  "mende.10species.10K":
    fq1: "../files/illumina_10species.10K.1.fq.gz"
    fq2: "../files/illumina_10species.10K.2.fq.gz"
  "mende.400species.10K":
    fq1: "../files/illumina_400species.10K.1.fq.gz"
    fq2: "../files/illumina_400species.10K.2.fq.gz"

run:
  ganon:
    "latest":
      dbs: 
        "bac_rs_refgen": "../build_test/ganon/latest/bac_rs_refgen/"
      fixed_args:
      binning_args:
        "--rel-cutoff": [0.25, 0.8]
      profiling_args:

  kmcp:
    "latest":
      dbs: 
        "bac_rs_refgen": "../build_test/kmcp/latest/bac_rs_refgen/"
      fixed_args:
      binning_args:
      profiling_args:
```

Verify run with `--dry-run`:

`snakemake -s metabench/classify.smk --configfile config/classify_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/classify.smk --configfile config/classify_test.yaml --cores 8 --use-conda`

If everything finished correctly, the following files will be created (`tree -A classify_test`):

```
classify_test/
├── ganon
│   └── latest
│       ├── mende.10species.10K
│       │   └── bac_rs_refgen
│       │       ├── default
│       │       │   ├── --rel-cutoff=0.25
│       │       │   │   ├── default.profiling.bench.json
│       │       │   │   ├── default.profiling.bench.tsv
│       │       │   │   ├── default.profiling.bioboxes.gz
│       │       │   │   └── default.profiling.log
│       │       │   ├── --rel-cutoff=0.25.binning.bench.json
│       │       │   ├── --rel-cutoff=0.25.binning.bench.tsv
│       │       │   ├── --rel-cutoff=0.25.binning.bioboxes.gz
│       │       │   ├── --rel-cutoff=0.25.binning.log
│       │       │   ├── --rel-cutoff=0.25.rep
│       │       │   ├── --rel-cutoff=0.8
│       │       │   │   ├── default.profiling.bench.json
│       │       │   │   ├── default.profiling.bench.tsv
│       │       │   │   ├── default.profiling.bioboxes.gz
│       │       │   │   └── default.profiling.log
│       │       │   ├── --rel-cutoff=0.8.binning.bench.json
│       │       │   ├── --rel-cutoff=0.8.binning.bench.tsv
│       │       │   ├── --rel-cutoff=0.8.binning.bioboxes.gz
│       │       │   ├── --rel-cutoff=0.8.binning.log
│       │       │   └── --rel-cutoff=0.8.rep
│       │       └── --max-fp=0.0001
│       │           ├── --rel-cutoff=0.25
│       │           │   ├── default.profiling.bench.json
│       │           │   ├── default.profiling.bench.tsv
│       │           │   ├── default.profiling.bioboxes.gz
│       │           │   └── default.profiling.log
│       │           ├── --rel-cutoff=0.25.binning.bench.json
│       │           ├── --rel-cutoff=0.25.binning.bench.tsv
│       │           ├── --rel-cutoff=0.25.binning.bioboxes.gz
│       │           ├── --rel-cutoff=0.25.binning.log
│       │           ├── --rel-cutoff=0.25.rep
│       │           ├── --rel-cutoff=0.8
│       │           │   ├── default.profiling.bench.json
│       │           │   ├── default.profiling.bench.tsv
│       │           │   ├── default.profiling.bioboxes.gz
│       │           │   └── default.profiling.log
│       │           ├── --rel-cutoff=0.8.binning.bench.json
│       │           ├── --rel-cutoff=0.8.binning.bench.tsv
│       │           ├── --rel-cutoff=0.8.binning.bioboxes.gz
│       │           ├── --rel-cutoff=0.8.binning.log
│       │           └── --rel-cutoff=0.8.rep
│       └── mende.400species.10K
│           └── bac_rs_refgen
│               ├── default
│               │   ├── --rel-cutoff=0.25
│               │   │   ├── default.profiling.bench.json
│               │   │   ├── default.profiling.bench.tsv
│               │   │   ├── default.profiling.bioboxes.gz
│               │   │   └── default.profiling.log
│               │   ├── --rel-cutoff=0.25.binning.bench.json
│               │   ├── --rel-cutoff=0.25.binning.bench.tsv
│               │   ├── --rel-cutoff=0.25.binning.bioboxes.gz
│               │   ├── --rel-cutoff=0.25.binning.log
│               │   ├── --rel-cutoff=0.25.rep
│               │   ├── --rel-cutoff=0.8
│               │   │   ├── default.profiling.bench.json
│               │   │   ├── default.profiling.bench.tsv
│               │   │   ├── default.profiling.bioboxes.gz
│               │   │   └── default.profiling.log
│               │   ├── --rel-cutoff=0.8.binning.bench.json
│               │   ├── --rel-cutoff=0.8.binning.bench.tsv
│               │   ├── --rel-cutoff=0.8.binning.bioboxes.gz
│               │   ├── --rel-cutoff=0.8.binning.log
│               │   └── --rel-cutoff=0.8.rep
│               └── --max-fp=0.0001
│                   ├── --rel-cutoff=0.25
│                   │   ├── default.profiling.bench.json
│                   │   ├── default.profiling.bench.tsv
│                   │   ├── default.profiling.bioboxes.gz
│                   │   └── default.profiling.log
│                   ├── --rel-cutoff=0.25.binning.bench.json
│                   ├── --rel-cutoff=0.25.binning.bench.tsv
│                   ├── --rel-cutoff=0.25.binning.bioboxes.gz
│                   ├── --rel-cutoff=0.25.binning.log
│                   ├── --rel-cutoff=0.25.rep
│                   ├── --rel-cutoff=0.8
│                   │   ├── default.profiling.bench.json
│                   │   ├── default.profiling.bench.tsv
│                   │   ├── default.profiling.bioboxes.gz
│                   │   └── default.profiling.log
│                   ├── --rel-cutoff=0.8.binning.bench.json
│                   ├── --rel-cutoff=0.8.binning.bench.tsv
│                   ├── --rel-cutoff=0.8.binning.bioboxes.gz
│                   ├── --rel-cutoff=0.8.binning.log
│                   └── --rel-cutoff=0.8.rep
└── kmcp
    └── latest
        ├── mende.10species.10K
        │   └── bac_rs_refgen
        │       └── default
        │           ├── default
        │           │   ├── default.profiling.bench.json
        │           │   ├── default.profiling.bench.tsv
        │           │   ├── default.profiling.bioboxes.gz
        │           │   └── default.profiling.log
        │           ├── default.binning.bench.json
        │           ├── default.binning.bench.tsv
        │           ├── default.binning.bioboxes.gz
        │           └── default.binning.log
        └── mende.400species.10K
            └── bac_rs_refgen
                └── default
                    ├── default
                    │   ├── default.profiling.bench.json
                    │   ├── default.profiling.bench.tsv
                    │   ├── default.profiling.bioboxes.gz
                    │   └── default.profiling.log
                    ├── default.binning.bench.json
                    ├── default.binning.bench.tsv
                    ├── default.binning.bioboxes.gz
                    └── default.binning.log
```

### Evaluations


```sh
mamba install snakemake
```