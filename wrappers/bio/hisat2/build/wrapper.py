__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

import re

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

extra = snakemake.params.get("extra", "")

fastas = snakemake.input.get("fastas") or snakemake.input.get("fasta")
seqs = snakemake.params.get("seqs")
assert (
    fastas is not None or seqs is not None
), "input/params: input fasta/fastas or params seqs is required"

if fastas:
    ref = [fastas] if isinstance(fastas, str) else fastas
    ref = ",".join([f"'{f}'" for f in ref])
    ref = f"-f {ref}"
else:
    ref = [seqs] if isinstance(seqs, str) else seqs
    ref = f"-c {ref}"

makedirs(snakemake.output)

shellcmd = (
    f"hisat2-build"
    f" -p {snakemake.threads}"
    f" {extra}"
    f" {ref}"
    f" '{snakemake.params.prefix}'"
    f" {log}"
)
shellcmd = re.sub(r"\s+", " ", shellcmd)
with open(snakemake.log[0], "wt") as log_fh:
    log_fh.write(f"{shellcmd}\n")
shell(shellcmd)
