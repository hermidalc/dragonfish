__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

db = snakemake.input.get("db")
assert db is not None, "input: db is a required input parameter"

fq = snakemake.input.get("fq")
assert fq is not None, "input: fq is a required input parameter"

extra = snakemake.params.get("extra", "")

shell(
    "diamond blastx"
    " --threads {snakemake.threads}"
    " --db {db}"
    " --query {fq}"
    " --out {snakemake.output[0]}"
    " {extra}"
    " {log}"
)
