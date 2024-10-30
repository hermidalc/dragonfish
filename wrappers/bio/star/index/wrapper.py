__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fastas = snakemake.input.get("fastas")
assert fastas is not None, "input: fastas is a required input parameter"

fastas = [fastas] if isinstance(fastas, str) else fastas
fastas = [f"'{f}'" for f in fastas]

gtf = snakemake.input.get("gtf", "")
if gtf:
    assert gtf.endswith((".gtf", ".gff3")), "input: gtf extension not .gtf/gff3"
    gtf = f"--sjdbGTFfile {gtf}"
    if gtf.endswith((".gff", ".gff3")):
        gtf += " --sjdbGTFtagExonParentTranscript Parent"

extra = snakemake.params.get("extra", "")

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --runMode genomeGenerate"
        " --genomeFastaFiles {fastas}"
        " --genomeDir {snakemake.output}"
        " --outTmpDir {tmp_dir}/tmp"
        " {gtf}"
        " {extra}"
        " {log}"
    )
