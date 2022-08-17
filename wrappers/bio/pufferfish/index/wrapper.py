__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

pufferfish = snakemake.params.get("pufferfish")
assert pufferfish is not None, "params: pufferfish is a required parameter"

ref = snakemake.input.get("ref")
assert ref is not None, "input: ref is a required input parameter"

out_dir = snakemake.output.get("out_dir", snakemake.output[0])

extra = snakemake.params.get("extra", "")

makedirs(out_dir)

tmp_base_dir = snakemake.params.get("tmp_dir", gettempdir())

with TemporaryDirectory(dir=tmp_base_dir) as tmp_dir:
    shell(
        "{pufferfish} index"
        " --threads {snakemake.threads}"
        " --ref {ref}"
        " --output {out_dir}"
        " --tmpdir {tmp_dir}"
        " {extra}"
        " {log}"
    )
