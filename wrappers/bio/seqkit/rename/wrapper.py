__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

fastas = snakemake.input.get("list_file")
if fastas is not None:
    fastas = f"--infile-list {fastas}"
else:
    fastas = snakemake.input

shell(
    "seqkit rename"
    " --threads {snakemake.threads}"
    " {extra}"
    " {fastas}"
    " --out-file {snakemake.output[0]}"
    " {log}"
)
