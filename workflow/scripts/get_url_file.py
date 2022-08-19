from hashlib import md5
from os.path import dirname
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlcleanup, urlretrieve

urlretrieve(snakemake.params[0], filename=snakemake.output[0])
urlcleanup()

try:
    md5_url = f"{snakemake.params[0]}.md5"
    md5_file = f"{snakemake.output[0]}.md5"
    urlretrieve(md5_url, filename=md5_file)
    urlcleanup()
except HTTPError:
    with open(snakemake.log[0], "w") as fh:
        fh.write(f"No {md5_url}\n")
else:
    with open(snakemake.log[0], "w") as fh:
        fh.write(f"Checking {md5_file}\n")
    with open(md5_file, "r") as fh:
        actual_md5, _ = fh.readline().strip().split(maxsplit=2)
    file_md5 = md5(Path(snakemake.output[0]).read_bytes()).hexdigest()
    assert (
        file_md5 == actual_md5
    ), f"File md5 {file_md5} doesn't match actual {actual_md5}"
