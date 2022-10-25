__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bam = snakemake.input.get("bam")
bam = snakemake.input[0] if bam is None else bam

flags = ""
read_count = snakemake.params.get("read_count", False)
if read_count:
    flags = "--count {snakemake.output[0]}"
extra = snakemake.params.get("extra")
if extra is not None:
    flags += " {extra}"

shell("seqkit bam --threads {snakemake.threads} {flags} {bam} {log}")
