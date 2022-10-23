__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

infiles = snakemake.input.get("list_file")
if infiles is not None:
    infiles = f"--infile-list {infiles}"
else:
    infiles = snakemake.input

flags = ""
id_regexp = snakemake.params.get("id_regexp")
if id_regexp is not None:
    flags += f"--id-regexp '{id_regexp}'"
extra = snakemake.params.get("extra")
if extra is not None:
    flags += f" {extra}"

shell(
    "seqkit seq"
    " --threads {snakemake.threads}"
    " {flags}"
    " {infiles}"
    " --out-file {snakemake.output[0]}"
    " {log}"
)
