from argparse import ArgumentParser
from os.path import basename
from urllib.parse import urljoin
from urllib.request import urlcleanup, urlretrieve

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--reports-url",
        type=str,
        default="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS",
        help="NCBI genomics assembly reports directory URL",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default="data/assembly_summary_genbank.txt",
        help="NCBI Genbank assembly summary file",
    )
    args = parser.parse_args()

    file_url = urljoin(f"{args.reports_url}/", basename(args.out_file))
    urlretrieve(file_url, filename=args.out_file)
    urlcleanup()
