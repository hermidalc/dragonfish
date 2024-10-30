__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

bam = snakemake.input.get("bam")
bam = snakemake.input[0] if bam is None else bam

gzip = snakemake.output[0].endswith(".gz")

flags = ""
read_count = snakemake.params.get("read_count", False)
if read_count:
    flags += "--count " + ("-" if gzip else snakemake.output[0]) + " --quiet-mode"

extra = snakemake.params.get("extra")
if extra:
    flags += f" {extra}"

pigz = f"| pigz -p {snakemake.threads} -c 1> {snakemake.output[0]}" if gzip else ""

log = snakemake.log_fmt_shell(stdout=False if gzip else True, stderr=True)

shell("seqkit bam --threads {snakemake.threads} {flags} {bam} {pigz} {log}")
