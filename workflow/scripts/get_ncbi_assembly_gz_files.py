from argparse import ArgumentParser
from os import makedirs
from os.path import join
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlcleanup, urlretrieve

import pandas as pd
from joblib import delayed, Parallel


def download_file(url, file, debug):
    if debug:
        print(url, flush=True)
    try:
        urlretrieve(url, filename=file)
        urlcleanup()
    except URLError as e:
        print(f"Missing {file}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--proteome-file",
        type=str,
        default="data/uniprot/uniprot_proteome_metadata.tsv",
        help="UniProt proteomes metadata file",
    )
    parser.add_argument(
        "--summary-file",
        type=str,
        default="data/genomes/assembly_summary_merged.txt",
        help="NCBI assembly summary merged file",
    )
    parser.add_argument(
        "--file-exts",
        nargs="+",
        type=str,
        default=["cds_from_genome.fna.gz", "genome.fna.gz", "genome.gtf.gz"],
        help="NCBI assembly file extensions",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="data/assembly",
        help="NCBI assembly output directory",
    )
    parser.add_argument("--n-jobs", type=int, default=-1, help="num parallel jobs")
    parser.add_argument(
        "--parallel-backend", type=str, default="loky", help="joblib parallel backend"
    )
    parser.add_argument("--verbose", type=int, default=1, help="Parallel verbosity")
    parser.add_argument("--debug", default=True, action="store_true", help="Debug")
    args = parser.parse_args()

    proteome_df = pd.read_csv(args.proteome_file, sep="\t", index_col="Proteome Id")

    summary_df = pd.read_csv(
        args.summary_file,
        sep="\t",
        index_col="assembly_accession",
        keep_default_na=False,
    )

    file_urls, files = [], []
    for acc in proteome_df["Genome assembly ID"].dropna():
        if acc in summary_df.index:
            dir_url = summary_df.loc[acc]["ftp_path"]
            dir_name = Path(dir_url).name
            for file_ext in args.file_exts:
                file_name = "_".join([dir_name, file_ext])
                file_urls.append(join(dir_url, file_name))
                files.append(join(args.out_dir, file_name))
        else:
            print(f"Missing {acc}", flush=True)

    makedirs(args.out_dir, mode=0o755, exist_ok=True)

    Parallel(n_jobs=args.n_jobs, backend=args.parallel_backend, verbose=args.verbose)(
        delayed(download_file)(url, file, args.debug)
        for url, file in zip(file_urls, files)
    )
