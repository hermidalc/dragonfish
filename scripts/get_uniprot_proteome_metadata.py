from argparse import ArgumentParser

import pandas as pd
from joblib import dump

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--api-url", type=str, help="Uniprot REST API URL")
    parser.add_argument(
        "--meta-file",
        type=str,
        default="data/uniprot_proteome_metadata.pkl",
        help="Dataframe out file",
    )
    args = parser.parse_args()

    metadata_df = pd.read_csv(args.api_url, sep="\t")
    dump(metadata_df, args.meta_file)
