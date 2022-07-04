from argparse import ArgumentParser
from ftplib import FTP

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--ftp-url", type=str, help="Uniprot REST API URL")
    parser.add_argument("--out-file", type=str, help="Dataframe out file")
    args = parser.parse_args()

    df = pd.read_csv(args.api_url, sep="\t")
    dump(df, args.out_file)
