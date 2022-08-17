__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

fastas = snakemake.input.get("fastas")
seqs = snakemake.params.get("seqs")
assert fastas is not None or seqs is not None, "input fastas or params seqs is required"

if fastas is not None:
    fastas = [fastas] if isinstance(fastas, str) else fastas
    fastas = ",".join([f"'{f}'" for f in fastas])
    seqs = ""
else:
    seqs = [seqs] if isinstance(seqs, str) else seqs
    seqs = f"-c {seqs}"

out_dir = snakemake.params.get("out_dir", "")
if out_dir:
    makedirs(out_dir)

shell(
    "hisat2-build"
    " -p {snakemake.threads}"
    " {extra}"
    " {fastas}"
    " {seqs}"
    " {out_dir}"
    " {log}"
)
