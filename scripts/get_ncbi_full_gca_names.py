from argparse import ArgumentParser
from os.path import split
from urllib.parse import urlparse

from entrezpy.conduit import Conduit
from joblib import load, dump

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--meta-file",
        type=str,
        default="data/uniprot_proteome_metadata.pkl",
        help="Uniprot proteomes metadata file",
    )
    parser.add_argument(
        "--name-file",
        type=str,
        default="data/nci_full_gca_names.pkl",
        help="Full GCA name file",
    )
    args = parser.parse_args()

full_gca_names = []
metadata_df = load(args.meta_file)
conduit = Conduit("hermidalc@pitt.edu")
for acc in metadata_df["Genome assembly ID"].dropna():
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
        print(acc, flush=True)
        num_uniq_ftp_dir_urls = len(set(ftp_dir_urls))
        if num_uniq_ftp_dir_urls > 1:
            raise ValueError(
                f"{acc} query returned {num_uniq_ftp_dir_urls:d} unique URLs"
            )
    ftp_dir_url = ftp_dir_urls[0]
    url_parts = urlparse(ftp_dir_url)
    full_gca_name = split(url_parts.path)[-1]
    full_gca_names.append(full_gca_name)

dump(full_gca_names, args.name_file)
