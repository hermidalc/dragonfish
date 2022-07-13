import gzip
from argparse import ArgumentParser
from os.path import splitext
from shutil import copyfileobj

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--file", type=str, help="Gzip file")
    args = parser.parse_args()

    with gzip.open(args.file, "rb") as f_in:
        with open(splitext(args.file)[0], "wb") as f_out:
            copyfileobj(f_in, f_out)
        f_in.close()
