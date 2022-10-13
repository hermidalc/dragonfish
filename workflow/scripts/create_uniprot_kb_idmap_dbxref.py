import vaex as vx

vx.settings.main.thread_count = snakemake.threads

idmap_df = vx.open(snakemake.input.idmap)
dbxref_df = vx.open(snakemake.input.dbxref)

idmap_dbxref_df = idmap_df.join(
    dbxref_df[dbxref_df.db == snakemake.params.db],
    how="inner",
    left_on="uniprot_id",
    right_on="uniprot_id",
    allow_duplication=True,
)
idmap_dbxref_df.drop("db", inplace=True)
idmap_dbxref_df.export_hdf5(snakemake.output[0])
