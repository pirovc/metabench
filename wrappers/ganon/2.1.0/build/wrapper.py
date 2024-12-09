"""Snakemake wrapper for ganon classify"""

__author__ = "Vitor C. Piro"
__copyright__ = "Copyright 2024, Vitor C. Piro"
__email__ = ""
__license__ = "MIT"

from snakemake.shell import shell
from os import path

# OUTPUT db prefix
db_prefix = path.splitext(snakemake.output[0])[0]

# PARAMS
input_extension = snakemake.params.get("input_extension", "")
if input_extension:
    input_extension = "--input-extension " + input_extension
taxonomy = snakemake.params.get("taxonomy", "")
if taxonomy:
    taxonomy = "--taxonomy " + taxonomy
taxonomy_files = snakemake.params.get("taxonomy_files", "")
if taxonomy_files:
    taxonomy_files = "--taxonomy-files " + taxonomy_files
ncbi_file_info = snakemake.params.get("ncbi_file_info", "")
if ncbi_file_info:
    ncbi_file_info = "--ncbi-file-info " + ncbi_file_info
extra = snakemake.params.get("extra", "")

# LOG
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "ganon build-custom "
    " --db-prefix {db_prefix}"
    " --input {snakemake.input}" 
    " --threads {snakemake.threads}"
    " {input_extension}"
    " {taxonomy}"
    " {taxonomy_files}"
    " {ncbi_file_info}"
    " {extra}"
    " {log}"
)
