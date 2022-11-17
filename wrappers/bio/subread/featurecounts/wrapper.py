__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from os import rename
from os.path import basename, exists, splitext

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

run_pigz = False
if snakemake.output[0].endswith(".gz"):
    output = splitext(snakemake.output[0])[0]
    run_pigz = True

shell(
    "featureCounts"
    " -T {snakemake.threads}"
    " -a {snakemake.input.gtf}"
    " {extra}"
    " -o {output}"
    " {snakemake.input.aligns}"
    " {log}"
)

if run_pigz:
    shell("pigz -p {snakemake.threads} {output}")

summary_file = f"{output}.summary"
if exists(summary_file):
    rename(summary_file, f"{basename(summary_file).partition('.')[0]}.summary")
