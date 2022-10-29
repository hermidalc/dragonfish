__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from os.path import splitext

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

cedar = snakemake.params.get("cedar")
assert cedar is not None, "params: cedar path is a required parameter"

assert snakemake.input[0].endswith(
    ("sam", "sam.gz", "pam")
), "input: file type required to be sam, sam.gz, or pam"

inflag = "--puffMapperOut" if snakemake.input[0].endswith(".pam") else "--sam"

run_pigz = False
if snakemake.output[0].endswith((".gz", ".GZ")):
    output = splitext(snakemake.output[0])[0]
    run_pigz = True
else:
    output = snakemake.output[0]

extra = snakemake.params.get("extra", "")

shell(
    "{cedar}"
    " --threads {snakemake.threads}"
    " {inflag} {snakemake.input[0]}"
    " --output {output}"
    " {extra}"
    " {log}"
)
if run_pigz:
    shell("pigz -p {snakemake.threads} {output}")
