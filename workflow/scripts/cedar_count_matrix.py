import pandas as pd


def load_cedar_df(cedar_file):
    cedar_df = pd.read_csv(
        cedar_file,
        sep="\t",
        index_col=0,
        usecols=[0, 3],
        engine="c",
        low_memory=False,
    )
    cedar_df.index.name = "tax_id"
    return cedar_df


merged_df = pd.DataFrame()
for cedar_file in snakemake.input:
    merged_df = pd.concat(
        [merged_df, load_cedar_df(cedar_file)], axis=0, verify_integrity=True
    )

merged_df.to_csv(snakemake.output[0], sep="\t")
