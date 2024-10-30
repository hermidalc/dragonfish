__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bam = snakemake.input.get("bam")
bam = snakemake.input[0] if bam is None else bam

extra = snakemake.params.get("extra", "")

qname = snakemake.input.get("qname")
if qname:
    extra = f"--qname-file {qname} {extra}"

shell(
    "samtools view"
    " --threads {snakemake.threads}"
    " {extra}"
    " --output {snakemake.output[0]}"
    " {bam}"
    " {log}"
)
