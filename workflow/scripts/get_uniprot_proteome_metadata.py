import pandas as pd

pd.read_csv(snakemake.input[0], sep="\t").to_csv(
    snakemake.output[0], sep="\t", index=False
)
