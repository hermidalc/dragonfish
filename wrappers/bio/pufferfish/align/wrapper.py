__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

pufferfish = snakemake.params.get("pufferfish")
assert pufferfish is not None, "params: pufferfish is a required parameter"

index = snakemake.input.get("index")
assert index is not None, "input: index is a required input parameter"

fq1 = snakemake.input.get("fq") or snakemake.input.get("fq1")
assert fq1 is not None, "input: fq/fq1 is a required input parameter"
fq2 = snakemake.input.get("fq2")
in_fqs = f"--mate1 {fq1} --mate2 {fq2}" if fq2 else f"--read {fq1}"

extra = snakemake.params.get("extra", "")

tmp_base_dir = snakemake.resources.get("tmpdir", gettempdir())

with TemporaryDirectory(dir=tmp_base_dir) as tmp_dir:
    shell(
        "{pufferfish} align"
        " --threads {snakemake.threads}"
        " --index {index}"
        " {in_fqs}"
        " --output {snakemake.output[0]}"
        " --tmpdir {tmp_dir}"
        " {extra}"
        " {log}"
    )
