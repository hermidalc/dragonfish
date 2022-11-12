from os.path import basename
from urllib.parse import urlparse

import pandas as pd

print("\nFiltering NCBI genome assemblies")

summary_df = pd.read_csv(
    snakemake.input.summary,
    sep="\t",
    index_col="assembly_accession",
    engine="c",
    low_memory=False,
)

proteome_df = pd.read_csv(
    snakemake.input.proteomes,
    sep="\t",
    index_col="Proteome Id",
    engine="c",
    low_memory=False,
)

accessions, genome_names = [], []
for acc in proteome_df["Genome assembly ID"]:
    if acc in summary_df.index:
        ftp_dir_url = summary_df.loc[acc]["ftp_path"]
        if pd.notna(ftp_dir_url):
            genome_name = basename(urlparse(ftp_dir_url).path)
            if genome_name in snakemake.params.skip:
                print(f"Skipping {genome_name}")
                continue
            accessions.append(acc)
            genome_names.append(genome_name)
        else:
            print(f"No FTP URL {acc}")
    else:
        print(f"Missing {acc}")

genome_name_series = pd.Series(genome_names)
dup_genome_names = genome_name_series[genome_name_series.duplicated()].values
assert (
    dup_genome_names.size == 0
), "Duplicate genome names found in UniProt Proteomes:\n{}".format(
    "\n".join(dup_genome_names)
)

summary_df = summary_df.loc[accessions]

print(f"\n{summary_df.shape[0]} filtered assemblies", flush=True)

summary_df.to_csv(snakemake.output[0], sep="\t")
