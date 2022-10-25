__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "diamond makedb"
    " --threads {snakemake.threads}"
    " --in {snakemake.input[0]}"
    " --db {snakemake.output[0]}"
    " {extra}"
    " {log}"
)