import pandas as pd

ref_df = pd.read_csv(snakemake.params.ref_url, sep="\t")
other_df = pd.read_csv(snakemake.params.other_url, sep="\t")

merged_df = pd.concat(
    [ref_df, other_df], axis=0, ignore_index=True, sort=False, verify_integrity=True
)
merged_df.dropna(subset="Genome assembly ID", inplace=True)
merged_df.reset_index(drop=True, inplace=True)

merged_df["Genome assembly ID NV"] = merged_df["Genome assembly ID"].str.split(
    ".", n=2, regex=False, expand=True
)[0]
merged_df.drop_duplicates(
    subset="Genome assembly ID NV", keep="first", inplace=True, ignore_index=True
)

merged_df.to_csv(snakemake.output[0], sep="\t", index=False)
