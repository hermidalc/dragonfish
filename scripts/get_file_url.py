from argparse import ArgumentParser
from urllib.request import urlcleanup, urlretrieve

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--file-url", type=str, help="File URL")
    parser.add_argument("--out-file", type=str, help="Output file")
    args = parser.parse_args()

    urlretrieve(args.file_url, filename=args.out_file)
    urlcleanup()
