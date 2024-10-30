__author__ = "Leandro C. Hermida"
__email__ = "leandro@leandrohermida.com"
__license__ = "MIT"

import re

from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


def get_format(path: str) -> str:
    """
    Return file format since HISAT2 reads files that
    could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    """
    if path.endswith((".gz", ".bz2")):
        return path.split(".")[-2].lower()
    return path.split(".")[-1].lower()


index = snakemake.params.get("idx")
assert index is not None, "params: idx is a required parameter"

n = len(snakemake.input.reads)
assert (
    n == 1 or n == 2
), "input: reads must have 1 (single-end) or 2 (paired-end) elements."

reads = snakemake.input.get("reads")
if isinstance(reads, str):
    reads = "-U {0}".format(reads)
elif len(reads) == 1:
    reads = "-U {0}".format(reads[0])
elif len(reads) == 2:
    reads = "-1 {0} -2 {1}".format(*reads)
else:
    raise RuntimeError(
        "Reads parameter must contain at least 1 and at most 2" " input files."
    )

extra = snakemake.params.get("extra", "")
if all(get_format(r) in ("fastq", "fq") for r in snakemake.input.reads):
    extra += " -q"
elif all(get_format(r) in ("fa", "mfa", "fasta") for r in snakemake.input.reads):
    extra += " -f"

unaligned = snakemake.output.get("unaligned")
if unaligned:
    un = (
        "un-gz"
        if unaligned.endswith((".gz", ".GZ"))
        else "un-bz2" if unaligned.endswith((".bz2", ".BZ2")) else "un"
    )
    extra += f" --{un} {unaligned}"

unpaired = snakemake.output.get("unpaired")
if unpaired:
    al = (
        "al-gz"
        if unpaired.endswith((".gz", ".GZ"))
        else "al-bz2" if unpaired.endswith((".bz2", ".BZ2")) else "al"
    )
    extra += f" --{al} {unpaired}"

unconcordant = snakemake.output.get("unconcordant")
if unconcordant:
    assert (
        isinstance(unconcordant, (list, tuple, set)) and len(unconcordant) == 2
    ), "input: unconcordant must be a list of two files"
    un_conc_groups_1 = re.split(
        r"(?:1|2)(\.(?:fq|fastq))(\.gz)?$", unconcordant[0], flags=re.IGNORECASE
    )
    un_conc_groups_2 = re.split(
        r"(?:1|2)(\.(?:fq|fastq))(\.gz)?$", unconcordant[1], flags=re.IGNORECASE
    )
    assert set(un_conc_groups_1) == set(
        un_conc_groups_2
    ), "input: unconcordant file name pair is invalid"
    un_conc_opt = (
        "un-conc-gz"
        if un_conc_groups_1[2].lower() == ".gz"
        else "un-conc-bz2" if un_conc_groups_1[2].lower() == ".bz2" else "un-conc"
    )
    un_conc_pat = f"{un_conc_groups_1[0]}%{un_conc_groups_1[1]}"
    if un_conc_groups_1[2].lower() in (".gz", ".bz2"):
        un_conc_pat = f"{un_conc_pat}{un_conc_groups_1[2]}"
    extra += f" --{un_conc_opt} {un_conc_pat}"

concordant = snakemake.output.get("concordant")
if concordant:
    assert (
        isinstance(concordant, (list, tuple, set)) and len(concordant) == 2
    ), "input: concordant must be a list of two files"
    al_conc_groups_1 = re.split(
        r"(?:1|2)(\.(?:fq|fastq))(\.gz)?$", concordant[0], flags=re.IGNORECASE
    )
    al_conc_groups_2 = re.split(
        r"(?:1|2)(\.(?:fq|fastq))(\.gz)?$", concordant[1], flags=re.IGNORECASE
    )
    assert set(al_conc_groups_1) == set(
        al_conc_groups_2
    ), "input: concordant file name pair is invalid"
    al_conc_opt = (
        "al-conc-gz"
        if al_conc_groups_1[2].lower() == ".gz"
        else "al-conc-bz2" if al_conc_groups_1[2].lower() == ".bz2" else "al-conc"
    )
    al_conc_pat = f"{al_conc_groups_1[0]}%{al_conc_groups_1[1]}"
    if al_conc_groups_1[2].lower() in (".gz", ".bz2"):
        al_conc_pat = f"{al_conc_pat}{al_conc_groups_1[2]}"
    extra += f" --{al_conc_opt} {al_conc_pat}"

hisat2_threads = snakemake.threads

# Setting parse_threads to false since samtools performs only
# bam compression. Thus the wrapper would use *twice* the amount
# of threads reserved by user otherwise.
samtools_opts = get_samtools_opts(snakemake, parse_threads=False)

output = snakemake.output.get("output")
if output is not None:
    assert output.lower().endswith((".bam", ".sam")), "output: must be SAM/BAM file"
    if output.lower().endswith(".bam"):
        hisat2_threads = hisat2_threads - 1
        if hisat2_threads < 1:
            raise ValueError(
                f"Wrapper expects at least two threads, got {hisat2_threads}"
            )
        output = f"| samtools view -h {samtools_opts}"
    else:
        output = f" -S {output}"
else:
    output = "> /dev/null"

log = snakemake.log_fmt_shell(
    stdout=False if output is None else True, stderr=True, append=True
)

shellcmd = (
    f"(hisat2"
    f" --threads {hisat2_threads}"
    f" {reads}"
    f" -x '{index}'"
    f" {extra}"
    f" {output}"
    f") {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
