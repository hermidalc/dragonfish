import re

import pandas as pd

print("\nGetting UniProt Proteomes metadata", flush=True)

ref_df = pd.read_csv(snakemake.params.ref_url, sep="\t", engine="c", low_memory=False)
other_df = pd.read_csv(
    snakemake.params.other_url, sep="\t", engine="c", low_memory=False
)

merged_df = pd.concat(
    [ref_df, other_df], axis=0, ignore_index=True, sort=False, verify_integrity=True
)
merged_df.dropna(subset="Genome assembly ID", inplace=True)
# merged_df.dropna(subset="Organism", inplace=True)

# keep newest assmbly version if there are duplicates
merged_df["Genome assembly ID NV"] = merged_df["Genome assembly ID"].str.split(
    ".", n=2, regex=False, expand=True
)[0]
merged_df = (
    merged_df.sort_values(by=["Genome assembly ID"])
    .drop_duplicates(subset="Genome assembly ID NV", keep="last", ignore_index=True)
    .sort_values(by=["Proteome Id"], ignore_index=True)
)

merged_df = merged_df.loc[
    (
        (
            merged_df["Taxonomic lineage"]
            .str.split("\s*,\s*", n=2, expand=True)[0]
            .str.capitalize()
            != "Eukaryota"
        )
        | (
            merged_df["Organism"]
            .str.capitalize()
            .str.startswith(tuple(snakemake.params.eukaryote_genera))
            .fillna(False)
        )
    )
    & ~merged_df["Organism"]
    .str.contains(snakemake.params.low_quality_pattern, regex=True, flags=re.IGNORECASE)
    .fillna(False)
]

print(f"\n{merged_df.shape[0]} UniProt proteomes", flush=True)

if snakemake.params.n_sample > 0:
    print(f"\nSampling {snakemake.params.n_sample} proteomes", flush=True)
    if snakemake.params.filter_domains:
        merged_df = merged_df.loc[
            ~merged_df["Taxonomic lineage"]
            .str.split("\s*,\s*", n=2, expand=True)[0]
            .str.capitalize()
            .isin(snakemake.params.filter_domains)
        ]
    merged_df = merged_df.sample(
        n=snakemake.params.n_sample,
        random_state=snakemake.params.random_seed,
        axis=0,
        ignore_index=True,
    )

merged_df.to_csv(snakemake.output[0], sep="\t", index=False)
