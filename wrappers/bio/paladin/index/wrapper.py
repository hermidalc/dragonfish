__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from os import remove
from os.path import basename, join
from shutil import copy

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

ref_type = snakemake.params.get("ref_type")
assert ref_type is not None, "input: ref_type is a required parameter"
assert int(ref_type) in range(1, 5), "input: ref_type must be 1..4"

extra = snakemake.params.get("extra", "")

gff = snakemake.input.get("gff", "")
if gff:
    assert ref_type == 1, "Error: ref_type must = 1 when providing gff file"

makedirs(snakemake.output[0])

ref = join(snakemake.output[0], basename(snakemake.input.ref).partition(".")[0])

copy(snakemake.input.ref, ref)

shell("paladin index -r{ref_type} {extra} {ref} {gff} {log}")

remove(ref)
