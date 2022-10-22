__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

min_orf_len = snakemake.params.get("min_orf_len", 250)

extra = snakemake.params.get("extra", "")

index = snakemake.input.get("index")
assert index is not None, "input: index is a required input parameter"

fq = snakemake.input.get("fq")
assert fq is not None, "input: fq is a required input parameter"

shell(
    "paladin align"
    " -t {snakemake.threads}"
    " -f {min_orf_len}"
    " {extra}"
    " {index}"
    " {fq}"
    " > {snakemake.output[0]}"
)
