__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from hashlib import md5
from os.path import join
from pathlib import Path
from urllib.request import urlretrieve, urlcleanup

regions2suffix = {
    "ALL": "chr_patch_hapl_scaff",
    "PRI": "primary_assembly",
}

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

protocol = snakemake.params.get("protocol")
assert protocol is not None, "params: protocol is a required parameter"
species = snakemake.params.get("species")
assert species is not None, "params: species is a required parameter"
species = species.lower()
release = snakemake.params.get("release")
assert release is not None, "params: release is a required  parameter"
regions = snakemake.params.get("regions")
assert regions is not None, "params: regions is a required parameter"
regions = regions.upper()
assert (
    regions == "CHR" or regions in regions2suffix
), f"params: regions {regions} is not valid"
annot_fmt = snakemake.params.get("annot_fmt")
assert annot_fmt is not None, "params: annot_fmt is a required parameter"
annot_fmt = annot_fmt.lower()
assert annot_fmt in ("gtf", "gff3"), f"params: annot_fmt {annot_fmt} is not valid"

base_url = f"{protocol}://ftp.ebi.ac.uk/pub/databases/gencode/"
release_url = join(base_url, f"Gencode_{species}", f"release_{release}")

gz_filename = (
    f"gencode.v{release}.annotation.{annot_fmt}.gz"
    if regions == "CHR"
    else f"gencode.v{release}.{regions2suffix[regions]}.annotation.{annot_fmt}.gz"
)
gz_file_url = join(release_url, gz_filename)

urlretrieve(gz_file_url, filename=snakemake.output[0])

name2md5 = {}
md5sums_url = join(release_url, "MD5SUMS")
md5sums_file = urlretrieve(md5sums_url)[0]
with open(md5sums_file, "r") as fh:
    for line in fh:
        md5sum, filename = line.strip().split(maxsplit=2)
        name2md5[filename] = md5sum

gz_file_md5 = md5(Path(snakemake.output[0]).read_bytes()).hexdigest()
assert (
    gz_file_md5 == name2md5[gz_filename]
), f"File md5sum {gz_file_md5} doesn't match actual {name2md5[gz_filename]}"

urlcleanup()
