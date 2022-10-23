__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

pufferfish = snakemake.params.get("pufferfish")
assert pufferfish is not None, "params: pufferfish path is a required parameter"

ref = snakemake.input.get("ref")
assert ref is not None, "input: ref FASTA is a required input parameter"

decoys = snakemake.input.get("decoys", "")
if decoys:
    decoys = "--decoys " + (decoys if isinstance(decoys, str) else f"<({decoys})")

extra = snakemake.params.get("extra", "")

tmp_base_dir = snakemake.resources.get("tmpdir", gettempdir())

with TemporaryDirectory(dir=tmp_base_dir) as tmp_dir:
    shell(
        "{pufferfish} index"
        " --threads {snakemake.threads}"
        " --ref {ref}"
        " --output {snakemake.output[0]}"
        " --tmpdir {tmp_dir}"
        " {decoys}"
        " {extra}"
        " {log}"
    )
