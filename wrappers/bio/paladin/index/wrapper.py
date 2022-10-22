__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from os import chdir

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

ref_type = snakemake.input.get("ref")
assert ref_type is not None, "input: ref_type is a required parameter"
assert int(ref_type) in range(1, 5), "input: ref_type must be 1..4"

extra = snakemake.params.get("extra", "")
annot = snakemake.input.get("gff", "")

makedirs(snakemake.output[0])
chdir(snakemake.output[0])

shell(
    "paladin index"
    " -r{ref_type}"
    " {extra}"
    " {snakemake.input[0]}"
    " {annot}"
    " {log}"
)
