from os.path import basename

import vaex as vx

file_basename = basename(snakemake.input[0]).partition(".")[0]

vx.from_csv(snakemake.input[0], convert=False, sep="\t", engine="c").export_hdf5(
    snakemake.output[0]
)
