from os.path import basename, splitext

import pandas as pd


def load_cedar_df(cedar_file):
    cedar_df = pd.read_csv(
        cedar_file,
        sep="\t",
        index_col=0,
        header=0,
        usecols=[0, 2],
        names=["tax_id", splitext(basename(cedar_file))[0]],
        engine="c",
        low_memory=False,
    )
    return cedar_df


merged_df = pd.DataFrame()
for cedar_file in snakemake.input:
    merged_df = pd.concat(
        [merged_df, load_cedar_df(cedar_file)], axis=1, verify_integrity=True
    )

merged_df.fillna(0, inplace=True)
merged_df.to_csv(snakemake.output[0], sep="\t")
