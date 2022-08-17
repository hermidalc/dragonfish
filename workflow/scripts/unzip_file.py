import gzip
from os.path import basename, dirname
from shutil import copyfileobj
from zipfile import ZipFile

try:
    if snakemake.input[0].endswith((".zip", ".ZIP")):
        with ZipFile(snakemake.input[0]) as zf:
            if len(snakemake.output) > 1:
                for member in snakemake.output:
                    zf.extract(basename(member), path=dirname(member))
            else:
                zf.extract(
                    basename(snakemake.output[0]), path=dirname(snakemake.output[0])
                )
    elif snakemake.input[0].endswith((".gz", ".GZ")):
        with gzip.open(snakemake.input[0], "rb") as f_in:
            with open(snakemake.output[0], "wb") as f_out:
                copyfileobj(f_in, f_out)
except Exception as e:
    print(f"Error {snakemake.input[0]}: {e}")
