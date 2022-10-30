import gzip

import pandas as pd

proteome_df = pd.read_csv(
    snakemake.input.proteomes,
    sep="\t",
    index_col="Proteome Id",
    engine="c",
    low_memory=False,
)
with open(snakemake.output[0], "wt") as out_fh:
    for i, gz_file in enumerate(sorted(snakemake.input.files, reverse=True)):
        with gzip.open(gz_file, "rt") as in_fh:
            if i == 0:
                out_fh.write(in_fh.readline())
            line = in_fh.readline()
            tax_id = line.split("\t")[2].strip()
            if tax_id in proteome_df["Organism Id"]:
                out_fh.write(line)
