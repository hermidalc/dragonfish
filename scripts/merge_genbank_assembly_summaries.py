from argparse import ArgumentParser
from shutil import copy
import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--summary-files",
        nargs="+",
        type=str,
        default=[
            "assembly_summary_genbank.txt",
            "assembly_summary_genbank_historical.txt",
        ],
        help="NCBI genomics assembly reports directory URL",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default="data/assembly_summary_merged.txt",
        help="NCBI Genbank assembly summary merged file",
    )
    args = parser.parse_args()

    if len(args.summary_files) > 1:
        merged_df = pd.DataFrame()
        for summary_file in args.summary_files:
            summary_df = pd.read_csv(
                summary_file, sep="\t", skiprows=1, index_col=0, na_values="na"
            )
            summary_df.index.name = "assembly_accession"
            merged_df = pd.concat([merged_df, summary_df], verify_integrity=True)
        merged_df.to_csv(args.out_file, sep="\t")
    else:
        copy(args.summary_files[0], args.out_file)
