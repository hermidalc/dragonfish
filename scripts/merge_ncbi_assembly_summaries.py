from argparse import ArgumentParser
from shutil import copy
import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--summary-files",
        nargs="+",
        type=str,
        help="NCBI assembly summary files",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default="data/assembly_summary_merged.txt",
        help="NCBI assembly summary merged file",
    )
    args = parser.parse_args()

    if len(args.summary_files) > 1:
        merged_df = pd.DataFrame()
        for summary_file in args.summary_files:
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
        merged_df.to_csv(args.out_file, sep="\t")
    else:
        copy(args.summary_files[0], args.out_file)
