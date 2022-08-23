import pandas as pd

print("\nGetting UniProt Proteomes metadata", flush=True)

ref_df = pd.read_csv(snakemake.params.ref_url, sep="\t")
other_df = pd.read_csv(snakemake.params.other_url, sep="\t")

merged_df = pd.concat(
    [ref_df, other_df], axis=0, ignore_index=True, sort=False, verify_integrity=True
)
merged_df.dropna(subset="Genome assembly ID", inplace=True)
# XXX: currently only one organism without a name
merged_df.dropna(subset="Organism", inplace=True)

merged_df["Genome assembly ID NV"] = merged_df["Genome assembly ID"].str.split(
    ".", n=2, regex=False, expand=True
)[0]
merged_df.drop_duplicates(
    subset="Genome assembly ID NV", keep="first", inplace=True, ignore_index=True
)

merged_df = merged_df[
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
    )
]

print(f"\n{merged_df.shape[0]} UniProt proteomes", flush=True)

if snakemake.params.n_sample > 0:
    print(f"\nSampling {snakemake.params.n_sample} proteomes", flush=True)
    if snakemake.params.bacteria_only:
        merged_df = merged_df[
            merged_df["Taxonomic lineage"]
            .str.split("\s*,\s*", n=2, expand=True)[0]
            .str.capitalize()
            == "Bacteria"
        ]
    merged_df = merged_df.sample(
        n=snakemake.params.n_sample,
        random_state=snakemake.params.random_seed,
        axis=0,
        ignore_index=True,
    )

merged_df.to_csv(snakemake.output[0], sep="\t", index=False)
