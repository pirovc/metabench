# MetaBench

MetaBench is a pipeline to continuously benchmark metagenomics analysis tools. It covers database construction (*build*), taxonomic *binning* and *profiling*.

It supports:
- multiple tools
	- in several release versions
- multiple databases
	- multiple range of parameters
- multiple samples
	- for binning and profiling both with multiple range of parameters
- multiple evaluation metrics and threshold
- performance benchmarks with optional repeats

It outputs:
- standardized JSON files for integration
- [bioboxes](https://github.com/bioboxes/rfc/tree/master/data-format) output files for binning and profiling
- interactive dashboard to analyze results

It requires:
- Reference sequences to build databases
- Set of reads in fastq and ground truth in the bioboxes format

Current configured tools:
- ganon v2.0.0
- kmcp v0.9.4
- metacache v2.3.1
- kraken2 v2.1.3
- bracken v2.9

## Installation and requirements

MetaBench is written in Snakemake and makes use of conda/mamba internally to install dependencies. It uses Bokeh to plot the interactive dashboard.

```sh
pip install randomname
mamba create -n metabench_env install snakemake pandas "bokeh==2.4.3"
source activate metabench_env
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
workdir: "example/build/"
threads: 8
repeat: 1

tools:
  ganon:
    "2.0.0": ""
  kmcp:
    "0.9.4": ""

dbs:
  "bac_rs_refgen":
    folder: "../../bac_rs/refgen/files/"
    extension: ".fna.gz"
    taxonomy: "ncbi"
    taxonomy_files: "../../bac_rs/refgen/taxdump.tar.gz"
    assembly_summary: "../../bac_rs/refgen/assembly_summary.txt"

run:
   ganon:
     "2.0.0":
       bac_rs_refgen:
         fixed_args:
           "--ncbi-file-info": "../../bac_rs/refgen/assembly_summary.txt"
         args:
           "--max-fp": [0.0001, ""]
   kmcp:
     "0.9.4":
       bac_rs_refgen:
         fixed_args:
         args:
```

In the example above MetaBench is set to build databases for 2 tools (ganon and kmcp). kmcp will run with default parameters only (no `args:`) and ganon will run with `--max-fp 0.0001` and default parameters.

Verify run with `--dry-run`:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda`

If everything finished correctly, the following files will be created (`tree -A example/build/`):

<pre>

```
example/build/
├── ganon
│   └── 2.0.0
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
    └── 0.9.4
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

</pre>

- `*.build.bench.json` contains the standardized metrics in JSON format. If `repeat > 1` in the config file, only the fastest run is selected.
- `*.build.bench.tsv` contains the raw benchmark metrics from Snakemake. If `repeat > 1` in the config file, one line for each run will be reported.
- `*.build.log` contains the STDOUT and STDERR from the run.
- `*.build.size.tsv` contains the size in bytes for the mandatory database files (`du --bytes`).

Obs: note that if no arguments are used in `args:` section of the configuration, the database folder/files will be named `default`. If parameters are used, databases are created based on them (`--max-fp 0.0001` -> `--max-fp=0.0001`, if more than one, connected by underscore "\_"). Any information provided in `fixed_args:` is not accounted for file/folder names.

Check the `config/build_example.yaml` for more examples on how to use the configuration file. Multiple databases, range of parameters and others can be configured to be executed in the same run.

### Classify (binning + profiling)

Classification includes both binning and profiling procedures. It requires databases (as created in the build process above) and one or more samples with single or paired `fastq` files.

Create `config/classify_test.yaml`:

```yaml
workdir: "example/classify/"
threads: 8
repeat: 1

tools:
  ganon:
    "2.0.0": ""
  kmcp:
    "0.9.4": ""

samples:
  "mende.10species.10K":
    fq1: "../../files/illumina_10species.10K.1.fq.gz"
    fq2: "../../files/illumina_10species.10K.2.fq.gz"

run:
  ganon:
    "2.0.0":
      dbs: 
        "bac_rs_refgen": "../../example/build/ganon/2.0.0/bac_rs_refgen/"
      fixed_args:
      binning_args:
        "--rel-cutoff": [0.25, 0.8]
      profiling_args:

  kmcp:
    "0.9.4":
      dbs: 
        "bac_rs_refgen": "../../example/build/kmcp/0.9.4/bac_rs_refgen/"
      fixed_args:
      binning_args:
      profiling_args:
```

Verify run with `--dry-run`:

`snakemake -s metabench/classify.smk --configfile config/classify_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/classify.smk --configfile config/classify_test.yaml --cores 8 --use-conda`

If everything finished correctly, the following files will be created (`tree -A example/classify/`):

<pre>

```
example/classify/
├── ganon
│   └── 2.0.0
│       └── mende.10species.10K
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
    └── 0.9.4
        └── mende.10species.10K
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

</pre>

- `*.profiling.bioboxes.gz` contains the standardized profiling output in bioboxes format.
- `*.binning.bioboxes.gz` contains the standardized binning output in bioboxes format.
- `*.bench.json` contains the runtime metrics in JSON format. If `repeat > 1` in the config file, only the fastest run is selected. This file will be updated in the next step (evaluation).
- `*.bench.tsv` contains the raw benchmark metrics from Snakemake. If `repeat > 1` in the config file, one line for each run will be reported.
- `*.log` contains the STDOUT and STDERR from the run.

### Evaluations

Evaluation will calculate metrics for binning and profiling procedures. It requires ground truth files for each sample

Create `config/evals_test.yaml`:

```yaml
workdir: "example/classify/"
threads: 8

samples:
  "mende.10species.10K":
    "binning": "../../files/illumina_10species.10K.binning.bioboxes.gz"
    "profiling": "../../files/illumina_10species.profile.bioboxes.gz"

# Optional, contents of the database for some metrics
dbs:
  "bac_rs_refgen": "../build/ganon/2.0.0/bac_rs_refgen/default/ganon_db.tax"

# Ranks to evaluate
ranks:
  - superkingdom
  - phylum
  - class
  - order
  - family
  - genus
  - species

taxonomy: "ncbi"
taxonomy_files: "../../bac_rs/refgen/taxdump.tar.gz"

# Set one or more thresholds for evaluation metrics [0-100]
threhsold_profiling:
  - 0

threhsold_binning:
  - 0
  - 0.05
  - 1
```

Verify run with `--dry-run`:

`snakemake -s metabench/evals.smk --configfile config/evals_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/evals.smk --configfile config/evals_test.yaml --cores 8 --use-conda`

If everything finished correctly, the following files will be created (`tree -A example/classify/`):

<pre>

```
example/classify/
├── ganon
│   └── 2.0.0
│       └── mende.10species.10K
│           └── bac_rs_refgen
│               ├── default
│               │   ├── --rel-cutoff=0.25
│               │   │   ├── default.profiling.bench.json
│               │   │   ├── default.profiling.bench.tsv
│               │   │   ├── default.profiling.bioboxes.gz
│               │   │   ├── default.profiling.evals.log
│               │   │   ├── default.profiling.log
│               │   │   └── default.profiling.updated_json
│               │   ├── --rel-cutoff=0.25.binning.bench.json
│               │   ├── --rel-cutoff=0.25.binning.bench.tsv
│               │   ├── --rel-cutoff=0.25.binning.bioboxes.gz
│               │   ├── --rel-cutoff=0.25.binning.evals.log
│               │   ├── --rel-cutoff=0.25.binning.log
│               │   ├── --rel-cutoff=0.25.binning.updated_json
│               │   ├── --rel-cutoff=0.25.rep
│               │   ├── --rel-cutoff=0.8
│               │   │   ├── default.profiling.bench.json
│               │   │   ├── default.profiling.bench.tsv
│               │   │   ├── default.profiling.bioboxes.gz
│               │   │   ├── default.profiling.evals.log
│               │   │   ├── default.profiling.log
│               │   │   └── default.profiling.updated_json
│               │   ├── --rel-cutoff=0.8.binning.bench.json
│               │   ├── --rel-cutoff=0.8.binning.bench.tsv
│               │   ├── --rel-cutoff=0.8.binning.bioboxes.gz
│               │   ├── --rel-cutoff=0.8.binning.evals.log
│               │   ├── --rel-cutoff=0.8.binning.log
│               │   ├── --rel-cutoff=0.8.binning.updated_json
│               │   └── --rel-cutoff=0.8.rep
│               └── --max-fp=0.0001
│                   ├── --rel-cutoff=0.25
│                   │   ├── default.profiling.bench.json
│                   │   ├── default.profiling.bench.tsv
│                   │   ├── default.profiling.bioboxes.gz
│                   │   ├── default.profiling.evals.log
│                   │   ├── default.profiling.log
│                   │   └── default.profiling.updated_json
│                   ├── --rel-cutoff=0.25.binning.bench.json
│                   ├── --rel-cutoff=0.25.binning.bench.tsv
│                   ├── --rel-cutoff=0.25.binning.bioboxes.gz
│                   ├── --rel-cutoff=0.25.binning.evals.log
│                   ├── --rel-cutoff=0.25.binning.log
│                   ├── --rel-cutoff=0.25.binning.updated_json
│                   ├── --rel-cutoff=0.25.rep
│                   ├── --rel-cutoff=0.8
│                   │   ├── default.profiling.bench.json
│                   │   ├── default.profiling.bench.tsv
│                   │   ├── default.profiling.bioboxes.gz
│                   │   ├── default.profiling.evals.log
│                   │   ├── default.profiling.log
│                   │   └── default.profiling.updated_json
│                   ├── --rel-cutoff=0.8.binning.bench.json
│                   ├── --rel-cutoff=0.8.binning.bench.tsv
│                   ├── --rel-cutoff=0.8.binning.bioboxes.gz
│                   ├── --rel-cutoff=0.8.binning.evals.log
│                   ├── --rel-cutoff=0.8.binning.log
│                   ├── --rel-cutoff=0.8.binning.updated_json
│                   └── --rel-cutoff=0.8.rep
└── kmcp
    └── 0.9.4
        └── mende.10species.10K
            └── bac_rs_refgen
                └── default
                    ├── default
                    │   ├── default.profiling.bench.json
                    │   ├── default.profiling.bench.tsv
                    │   ├── default.profiling.bioboxes.gz
                    │   ├── default.profiling.evals.log
                    │   ├── default.profiling.log
                    │   └── default.profiling.updated_json
                    ├── default.binning.bench.json
                    ├── default.binning.bench.tsv
                    ├── default.binning.bioboxes.gz
                    ├── default.binning.evals.log
                    ├── default.binning.log
                    └── default.binning.updated_json
```

</pre>

- `*.bench.json` were updated with evaluation metrics.
- `*.evals.log` contains the STDOUT and STDERR from the evaluation run.
- `*.updated_json` flag file to set the evaluation is over and updated the correspondent .json file

### Plotting

Finally, to visualize the benchmark, plot the results:

```sh
scripts/plot.py -i example/ --output example_dashboard.html
```

Open the `example_dashboard.html` in your browser and explore the results.
