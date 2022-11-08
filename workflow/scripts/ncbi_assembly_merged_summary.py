import pandas as pd


def load_summary_df(summary_file):
    summary_df = pd.read_csv(
        summary_file,
        sep="\t",
        skiprows=1,
        index_col=0,
        na_values=("na", "NA"),
        engine="c",
        low_memory=False,
    )
    summary_df.index.name = "assembly_accession"
    return summary_df


merged_df = pd.DataFrame()
for summary_file in snakemake.input:
    merged_df = pd.concat(
        [merged_df, load_summary_df(summary_file)],
        axis=0,
        ignore_index=True,
        verify_integrity=True,
    )

merged_df.to_csv(snakemake.output[0], sep="\t")
