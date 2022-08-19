__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

pattern = snakemake.params.get("pattern")
assert pattern is not None, "params: pattern is a required parameter"

replacement = snakemake.params.get("replacement")
if replacement:
    replacement = f"--replacement '{replacement}'"

extra = snakemake.params.get("extra", "")

shell(
    "seqkit replace"
    " --threads {snakemake.threads}"
    " --pattern {pattern:q}"
    " {replacement}"
    " {extra}"
    " {snakemake.input}"
    " --out-file {snakemake.output}"
    " {log}"
)
