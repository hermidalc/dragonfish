import zipfile
from argparse import ArgumentParser
from os.path import basename, dirname

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--file", type=str, help="Zip file archive")
    parser.add_argument(
        "--members", type=str, nargs="+", help="Archive members to unzip"
    )
    parser.add_argument("--out-dir", type=str, help="Output directory")
    args = parser.parse_args()

    out_dir = dirname(args.file) if not args.out_dir else args.out_dir

    zfile = zipfile.ZipFile(args.file)
    for member in args.members:
        zfile.extract(basename(member), path=out_dir)
    zfile.close()
