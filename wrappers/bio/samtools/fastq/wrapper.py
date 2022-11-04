__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

pigz = (
    f"| pigz -p {snakemake.threads} -c 1> {snakemake.output[0]}"
    if snakemake.output[0].endswith(".gz")
    else ""
)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "samtools fastq --threads {snakemake.threads} {extra} {snakemake.input[0]} {pigz} {log}"
)
