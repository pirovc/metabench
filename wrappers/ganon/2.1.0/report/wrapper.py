"""Snakemake wrapper for ganon report"""

__author__ = "Vitor C. Piro"
__copyright__ = "Copyright 2024, Vitor C. Piro"
__email__ = ""
__license__ = "MIT"

from snakemake.shell import shell
from os import path, rename

# OUTPUT
output_prefix = path.splitext(snakemake.output[0])[0]

print(output_prefix)

# PARAMS
ranks = snakemake.params.get("ranks", "")
if ranks:
    ranks = "--ranks " + ranks

extra = snakemake.params.get("extra", "")

# LOG
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "ganon report "
    " --input {snakemake.input.rep}"
    " --db-prefix {snakemake.input.db_tax}"
    " --output-prefix {output_prefix}"
    " --output-format bioboxes"
    " {ranks}"
    " {extra}"
    " {log}"
)

rename(output_prefix + ".tre", snakemake.output[0])
