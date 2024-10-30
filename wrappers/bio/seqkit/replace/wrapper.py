__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

infiles = snakemake.input.get("list_file")
if infiles:
    infiles = f"--infile-list {infiles}"
else:
    infiles = snakemake.input

flags = ""
pattern = snakemake.params.get("pattern")
if pattern:
    flags += f"--pattern '{pattern}'"
    replacement = snakemake.params.get("replacement")
    if replacement:
        flags += f" --replacement '{replacement}'"
id_regexp = snakemake.params.get("id_regexp")
if id_regexp:
    flags += f"--id-regexp '{id_regexp}'"
extra = snakemake.params.get("extra")
if extra:
    flags += f" {extra}"

shell(
    "seqkit replace"
    " --threads {snakemake.threads}"
    " {flags}"
    " {infiles}"
    " --out-file {snakemake.output[0]}"
    " {log}"
)
