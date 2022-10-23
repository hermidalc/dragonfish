__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

cedar = snakemake.params.get("cedar")
assert cedar is not None, "params: cedar path is a required parameter"

extra = snakemake.params.get("extra", "")

shell(
    "{cedar}"
    " --threads {snakemake.threads}"
    " --puffMapperOut {snakemake.input[0]}"
    " --output {snakemake.output[0]}"
    " {extra}"
    " {log}"
)
