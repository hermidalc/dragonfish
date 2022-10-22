__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

extra = snakemake.params.get("extra", "")

shell("samtools view {extra} {snakemake.input[0]} {log}")
