from hashlib import md5
from os.path import basename, dirname, join
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import urlcleanup, urlopen, urlretrieve

import pandas as pd
from joblib import delayed, Parallel
from snakemake.utils import makedirs

md5_file_name = "md5checksums.txt"


def download_file(url, file, debug):
    if debug:
        print(url, flush=True)
    makedirs(dirname(file))
    urlretrieve(url, filename=file)
    urlcleanup()


def check_md5(file):
    md5_file = join(dirname(file), md5_file_name)
    md5_df = pd.read_csv(
        md5_file, delim_whitespace=True, header=None, names=["md5", "name"]
    )
    md5_df["name"] = md5_df["name"].apply(lambda n: basename(n))
    actual_md5 = md5_df["md5"][md5_df["name"] == basename(file)].values[0]
    file_md5 = md5(Path(file).read_bytes()).hexdigest()
    assert (
        file_md5 == actual_md5
    ), f"File md5 {file_md5} doesn't match actual {actual_md5}"


proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

summary_df = pd.read_csv(
    snakemake.input.summary_file,
    sep="\t",
    index_col="assembly_accession",
    low_memory=False,
)

genome_names = set()
gz_file_urls, gz_files = [], []
md5_file_urls, md5_files = [], []
for acc in proteome_df["Genome assembly ID"].dropna():
    if acc in summary_df.index:
        ftp_dir_url = summary_df.loc[acc]["ftp_path"]
        if pd.notna(ftp_dir_url):
            genome_name = basename(urlparse(ftp_dir_url).path)
            assert (
                genome_name not in genome_names
            ), f"{genome_name} appears more than once in UniProt proteome metadata"
            genome_names.add(genome_name)
            file_group_exists = True
            for gz_file_ext in snakemake.params.file_exts:
                gz_file_name = "_".join([genome_name, gz_file_ext])
                gz_file_url = join(ftp_dir_url, gz_file_name)
                try:
                    urlopen(gz_file_url)
                except HTTPError:
                    file_group_exists = False
                    break
            if file_group_exists:
                genome_out_dir = join(snakemake.output[0], genome_name)
                for gz_file_ext in snakemake.params.file_exts:
                    gz_file_name = "_".join([genome_name, gz_file_ext])
                    gz_file_url = join(ftp_dir_url, gz_file_name)
                    gz_file = join(genome_out_dir, gz_file_name)
                    gz_file_urls.append(gz_file_url)
                    gz_files.append(gz_file)
                md5_file = join(genome_out_dir, md5_file_name)
                md5_file_urls.append(join(ftp_dir_url, md5_file_name))
                md5_files.append(md5_file)
            else:
                print(f"Missing {acc}", flush=True)
        else:
            print(f"Missing {acc}", flush=True)
    else:
        print(f"Missing {acc}", flush=True)

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.parallel_backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(url, file, snakemake.params.debug)
    for url, file in zip(gz_file_urls + md5_file_urls, gz_files + md5_files)
)

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.parallel_backend,
    verbose=snakemake.params.verbosity,
)(delayed(check_md5)(gz_file) for gz_file in gz_files)
