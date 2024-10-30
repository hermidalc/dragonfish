__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

import re
from os.path import join
from shutil import rmtree
from tempfile import gettempdir, TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq") or snakemake.input.get("fq1")
assert fq1 is not None, "input: fq/fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
fq1_str = ",".join(fq1)

fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
    )
    assert len(fq1) == len(fq2), "input: equal number of files required for fq1 and fq2"
    fq2_str = ",".join(fq2)
    fq_str = f"'{fq1_str}' '{fq2_str}'"
else:
    fq_str = f"'{fq1_str}'"

index = snakemake.input.get("index", snakemake.params.get("index", "GenomeDir/"))

out_dir = snakemake.params.get("out_dir")
assert out_dir is not None, "params: out_dir is a required parameter"

readcmd = "--readFilesCommand zcat" if fq1[0].endswith(".gz") else ""

gtf = snakemake.input.get("gtf", "")
if gtf:
    assert gtf.endswith((".gtf", ".gff3")), "input: gtf extension not .gtf/gff3"
    gtf = f"--sjdbGTFfile {gtf}"
    if gtf.endswith((".gff", ".gff3")):
        gtf += " --sjdbGTFtagExonParentTranscript Parent"

read_length = snakemake.params.get("read_length")
if read_length is None:
    read_length_file = snakemake.input.get("read_length")
    assert (
        read_length_file is not None
    ), "input/params: read_length is a required parameter"
    with open(read_length_file, "r") as fh:
        read_length = re.sub("\D+", "", fh.readline())

sjdb_overhang = int(read_length) - 1 if read_length else 100

extra = snakemake.params.get("extra", "")
sj = snakemake.input.get("sj")
if sj:
    extra += f" --sjdbFileChrStartEnd {sj}"

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    shell(
        "STAR"
        " --runThreadN {snakemake.threads}"
        " --readFilesIn {fq_str}"
        " --genomeDir {index}"
        " --outFileNamePrefix {out_dir}/"
        " --outTmpDir {tmp_dir}/tmp"
        " --sjdbOverhang {sjdb_overhang}"
        " {readcmd}"
        " {gtf}"
        " {extra}"
        " {log}"
    )

rmtree(join(out_dir, "_STARgenome"), ignore_errors=True)
