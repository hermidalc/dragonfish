__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

if snakemake.input[0].endswith(".gz"):
    shell(
        "pigz -dc {snakemake.input[0]} |"
        " gffread"
        " {extra}"
        " -o {snakemake.output[0]}"
        " {log}"
    )
else:
    shell(
        "gffread"
        " {extra}"
        " -o {snakemake.output[0]}"
        " {snakemake.input[0]}"
        " {log}"
    )
