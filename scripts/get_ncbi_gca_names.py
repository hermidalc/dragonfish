from argparse import ArgumentParser
from os.path import split
from urllib.parse import urlparse

import pandas as pd
from entrezpy.conduit import Conduit
from joblib import delayed, Parallel


def get_ncbi_gca_name(acc, debug):
    if debug:
        print(acc, flush=True)
    conduit = Conduit("hermidalc@pitt.edu")
    pipe = conduit.new_pipeline()
    pid = pipe.add_search(
        {"db": "assembly", "term": acc, "field": "Assembly Accession"}
    )
    pid = pipe.add_summary(dependency=pid)
    analyzer = conduit.run(pipe)
    ftp_dir_urls = [
        summary["ftppath_genbank"]
        for summary in analyzer.get_result().summaries.values()
    ]
    if len(ftp_dir_urls) > 1:
        num_uniq_ftp_dir_urls = len(set(ftp_dir_urls))
        if num_uniq_ftp_dir_urls > 1:
            raise ValueError(
                f"{acc} query returned {num_uniq_ftp_dir_urls:d} unique URLs"
            )
    ftp_dir_url = ftp_dir_urls[0]
    url_parts = urlparse(ftp_dir_url)
    gca_name = split(url_parts.path)[-1]
    return gca_name


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--meta-file",
        type=str,
        default="data/uniprot_proteome_metadata.tsv",
        help="UniProt proteomes metadata TSV file",
    )
    parser.add_argument(
        "--name-file",
        type=str,
        default="data/ncbi_gca_names.tsv",
        help="NCBI GCA name TSV file",
    )
    parser.add_argument("--n-jobs", type=int, default=-1, help="num parallel jobs")
    parser.add_argument(
        "--parallel-backend", type=str, default="loky", help="joblib parallel backend"
    )
    parser.add_argument("--verbose", type=int, default=0, help="Parallel verbosity")
    parser.add_argument("--debug", default=False, action="store_true", help="Debug")
    args = parser.parse_args()

    metadata_df = pd.read_csv(args.meta_file, sep="\t")

    gca_names = Parallel(
        n_jobs=args.n_jobs, backend=args.parallel_backend, verbose=args.verbose
    )(
        delayed(get_ncbi_gca_name)(acc, args.debug)
        for acc in metadata_df["Genome assembly ID"].dropna()
    )
    pd.DataFrame(gca_names).to_csv(args.name_file, sep="\t", index=False, header=False)
