import vaex as vx

df = vx.open_many(snakemake.input)
df.export_hdf5(
    snakemake.output[0], column_count=df.shape[1], writer_threads=df.shape[1]
)
