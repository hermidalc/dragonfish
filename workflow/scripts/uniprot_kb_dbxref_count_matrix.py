import os

import vaex as vx

os.environ["MKL_NUM_THREADS"] = str(snakemake.threads)
os.environ["NUMBA_NUM_THREADS"] = str(snakemake.threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(snakemake.threads)
os.environ["OPENBLAS_NUM_THREADS"] = str(snakemake.threads)

vx.settings.main.thread_count = snakemake.threads
vx.settings.main.thread_count_io = snakemake.threads

genbank_counts_df = vx.open(snakemake.input.counts)
genbank_counts_df.rename(genbank_counts_df.column_names[0], "genbank_id")

idmap_dbxref_df = vx.open(snakemake.input.idmap)

dbxref_counts_df = idmap_dbxref_df.join(
    genbank_counts_df,
    how="inner",
    left_on="genbank_id",
    right_on="genbank_id",
    allow_duplication=True,
)

dbxref_counts_df.drop(["genbank_id", "uniprot_id"], inplace=True)

dbxref_counts_df.groupby(by="db_id", agg="sum", sort=True)

dbxref_counts_df.export_csv_arrow(snakemake.output[0], chunk_size=int(1e7))
