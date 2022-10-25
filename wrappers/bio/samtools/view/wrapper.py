__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

bam = snakemake.input.get("bam")
bam = snakemake.input[0] if bam is None else bam

extra = snakemake.params.get("extra", "")
qname = snakemake.params.get("qname")
if qname is not None:
    extra = "--qname-file {qname} {extra}"

shell(
    "samtools view"
    " --threads {snakemake.threads}"
    " {extra}"
    " --output {snakemake.output[0]}"
    " {bam}"
    " {log}"
)
