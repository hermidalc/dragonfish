from os import makedirs
from os.path import dirname
from urllib.request import urlcleanup, urlretrieve

makedirs(dirname(snakemake.output[0]), mode=0o755, exist_ok=True)

urlretrieve(snakemake.input[0], filename=snakemake.output[0])
urlcleanup()
