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
workdir: "databases/"
threads: 8
repeat: 1

tools:
  ganon:
    "2.0.0": ""
  kmcp:
    "0.9.2": ""

dbs:
  "bac_rs_refgen":
    folder: "../bac_rs/refgen/files/"
    extension: ".fna.gz"
    taxonomy: "ncbi"
    taxonomy_files: "../bac_rs/refgen/taxdump.tar.gz"
    assembly_summary: "../bac_rs/refgen/assembly_summary.txt"

run:
   ganon:
     "2.0.0":
       bac_rs_refgen:
         fixed_args:
           "--ncbi-file-info": "../bac_rs/refgen/assembly_summary.txt"
         args:
   kmcp:
     "0.9.2":
       bac_rs_refgen:
         fixed_args:
         args:
```

In the example above MetaBench is set to build databases for 2 tools (ganon and kmco) with default parameters for the previously downloaded files:

Verify run with `--dry-run`:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda --dry-run`

Run it:

`snakemake -s metabench/build.smk --configfile config/build_test.yaml --cores 8 --use-conda`

If everything ran correctly, the database files will look something like:

```
$ tree -A databases

	databases/
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
	│           └── default.build.size.tsv
	└── kmcp
	    └── 0.9.2
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

	11 directories, 27 files
```

Check the `config/build_example.yaml` for more examples on how to use the configuration file. Multiple databases or configuration can be executed in the same run.

### Classify (binning + profiling)


### Evaluations


```sh
mamba install snakemake
```