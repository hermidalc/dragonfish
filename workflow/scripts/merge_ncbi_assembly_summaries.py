from shutil import copy

import pandas as pd

if len(args.summary_files) > 1:
    merged_df = pd.DataFrame()
    for summary_file in snakemake.input[0]:
        summary_df = pd.read_csv(
            summary_file,
            sep="\t",
            skiprows=1,
            index_col=0,
            keep_default_na=False,
        )
        summary_df.index.name = "assembly_accession"
        summary_df.fillna("na", inplace=True)
        summary_df.fillna("NA", inplace=True)
        merged_df = pd.concat([merged_df, summary_df], verify_integrity=True)
    merged_df.to_csv(snakemake.output[0], sep="\t")
else:
    copy(snakemake.input[0], snakemake.output[0])
