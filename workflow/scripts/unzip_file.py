import gzip
from os.path import basename, dirname, splitext
from shutil import copyfileobj
from zipfile import ZipFile

if snakemake.input[0].endswith((".zip", ".ZIP")):
    with ZipFile(snakemake.input[0]) as zf:
        for member in snakemake.output[0]:
            zf.extract(basename(member), path=dirname(member))
elif snakemake.input[0].endswith((".gz", ".GZ")):
    with gzip.open(snakemake.input[0], "rb") as f_in:
        with open(splitext(snakemake.output[0])[0], "wb") as f_out:
            copyfileobj(f_in, f_out)
