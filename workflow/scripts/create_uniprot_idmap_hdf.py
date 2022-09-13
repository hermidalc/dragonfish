from os.path import basename, join

import vaex as vx
from snakemake.utils import makedirs

makedirs(snakemake.output[0])

file_basename = basename(snakemake.input[0]).partition(".")[0]

for n, df in enumerate(
    vx.from_csv(
        snakemake.input[0],
        chunk_size=float(snakemake.params.chunk_size),
        sep="\t",
        names=["uniparc", "genbank"],
        usecols=[10, 17],
        engine="c",
    ),
    start=1,
):
    df.export_hdf5(join(snakemake.output[0], f"{file_basename}_{n:02}.hdf5"))
