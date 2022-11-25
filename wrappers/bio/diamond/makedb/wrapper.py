__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from os.path import splitext

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fastas = snakemake.input.get("fastas")
assert fastas is not None, "input: fastas is a required input parameter"

extra = snakemake.params.get("extra", "")

taxonmap = snakemake.input.get("taxonmap")
if taxonmap:
    extra += f" --taxonmap {taxonmap}"
taxonnodes = snakemake.input.get("taxonnodes")
if taxonnodes:
    extra += f" --taxonnodes {taxonnodes}"
taxonnames = snakemake.input.get("taxonnames")
if taxonnames:
    extra += f" --taxonnames {taxonnames}"

output = splitext(snakemake.output[0])[0]

if isinstance(fastas, str):
    shell(
        "diamond makedb"
        " --threads {snakemake.threads}"
        " --in {fastas}"
        " --db {output}"
        " {extra}"
        " {log}"
    )
else:
    shell(
        "pigz -dc {fastas} | diamond makedb"
        " --threads {snakemake.threads}"
        " --db {output}"
        " {extra}"
        " {log}"
    )
