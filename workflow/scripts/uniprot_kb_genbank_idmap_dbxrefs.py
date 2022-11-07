import os

import vaex as vx

os.environ["MKL_NUM_THREADS"] = str(snakemake.threads)
os.environ["NUMBA_NUM_THREADS"] = str(snakemake.threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(snakemake.threads)
os.environ["OPENBLAS_NUM_THREADS"] = str(snakemake.threads)

vx.settings.main.thread_count = snakemake.threads
vx.settings.main.thread_count_io = snakemake.threads

idmap_df = vx.open(snakemake.input.idmap)
dbxref_df = vx.open(snakemake.input.dbxref)

idmap_dbxref_df = idmap_df.join(
    dbxref_df[dbxref_df.db == snakemake.params.db].extract(),
    how="inner",
    left_on="uniprot_id",
    right_on="uniprot_id",
    allow_duplication=True,
)
idmap_dbxref_df.drop("db", inplace=True)
idmap_dbxref_df.export_hdf5(snakemake.output[0], column_count=3, writer_threads=3)
