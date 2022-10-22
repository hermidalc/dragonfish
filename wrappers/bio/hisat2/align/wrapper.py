__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from os.path import normpath, sep

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

index = snakemake.input.get("index")
assert index is not None, "input: index is a required input parameter"

fq1 = snakemake.input.get("fq") or snakemake.input.get("fq1")
assert fq1 is not None, "input: fq/fq1 is a required input parameter"

fq2 = snakemake.input.get("fq2")
fqs = f"-1 '{fq1}' -2 '{fq}'" if fq2 else f"-U '{fq1}'"

extra = snakemake.params.get("extra", "")

shell("hisat2 -p {snakemake.threads} -x {index} {fqs} {extra} {log}")
