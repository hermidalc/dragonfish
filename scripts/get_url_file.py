from argparse import ArgumentParser
from os import makedirs
from os.path import dirname
from urllib.request import urlcleanup, urlretrieve

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--file-url", type=str, help="File URL")
    parser.add_argument("--out-file", type=str, help="Output file")
    args = parser.parse_args()

    out_dir = dirname(args.out_file)
    makedirs(out_dir, mode=0o755, exist_ok=True)

    urlretrieve(args.file_url, filename=args.out_file)
    urlcleanup()
