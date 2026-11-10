"""Snakemake wrapper for ganon classify"""

__author__ = "Vitor C. Piro"
__copyright__ = "Copyright 2026, Vitor C. Piro"
__email__ = ""
__license__ = "MIT"

from snakemake.shell import shell
from os import path, getcwd

# INPUT reads
single_reads = snakemake.input.get("single_reads", "")
if single_reads:
    single_reads = "--single-reads " + (" ".join(single_reads) if isinstance(single_reads, list) else single_reads)
paired_reads = snakemake.input.get("paired_reads", "")
if paired_reads:
    paired_reads = "--paired-reads " + (" ".join(paired_reads) if isinstance(paired_reads, list) else paired_reads)

batch_reads = snakemake.input.get("batch_reads", "")
if batch_reads:
    batch_reads = f"--batch-reads {batch_reads}"
    single_reads = ""
    paired_reads = ""

assert single_reads + paired_reads + batch_reads, "input-> single_reads or paired_reads required"

# INPUT db
dbs = snakemake.input.get("dbs")
db_prefixes = []
for db in dbs:
    assert db.endswith("ibf"), "input-> db should be an .ibf or .hibf file"
    db_prefixes.append(path.splitext(db)[0])

# OUTPUT
if batch_reads:
    output_prefix = getcwd() + "/"
else:
    out = {}
    ext = ["rep", "one", "all", "tre", "unc"]
    for e in ext:
        out[e] = snakemake.output.get(e)
        if out[e]:
            output_prefix = path.splitext(out[e])[0]

# PARAMS
extra = snakemake.params.get("extra", "")

# LOG
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "ganon classify "
    " --db-prefix {db_prefixes}"
    " --output-prefix {output_prefix}" 
    " --threads {snakemake.threads}"
    " --output-one"
    " {single_reads}"
    " {paired_reads}"
    " {batch_reads}"
    " {extra}"
    " {log}"
)
