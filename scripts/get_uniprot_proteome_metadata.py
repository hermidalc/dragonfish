from argparse import ArgumentParser

import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--api-url", type=str, help="Uniprot REST API URL")
    parser.add_argument(
        "--out-file",
        type=str,
        default="data/uniprot_proteome_metadata.tsv",
        help="UniProt proteomes metadata TSV file",
    )
    args = parser.parse_args()

    metadata_df = pd.read_csv(args.api_url, sep="\t")
    metadata_df.to_csv(args.out_file, sep="\t", index=False)
