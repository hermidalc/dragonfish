import vaex as vx

vx.open_many(snakemake.input).export_hdf5(snakemake.output[0])