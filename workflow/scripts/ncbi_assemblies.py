from hashlib import md5
from os import listdir, remove
from os.path import basename, dirname, exists, join
from pathlib import Path
from shutil import rmtree
from time import sleep
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import urlcleanup, urlretrieve

import pandas as pd
from joblib import delayed, Parallel
from snakemake.utils import makedirs


def download_file(url, file, retries, retry_wait):
    makedirs(dirname(file))
    while retries > 0:
        print(f"Downloading {url}", flush=True)
        try:
            urlretrieve(url, filename=file)
        except HTTPError as e:
            if e.code == 404:
                print(f"Skipped: {url}: {e}", flush=True)
                break
            if exists(file):
                remove(file)
            print(f"Retrying: {url}: {e}", flush=True)
            retries -= 1
            sleep(retry_wait)
        except Exception as e:
            if exists(file):
                remove(file)
            print(f"Retrying: {url}: {e}", flush=True)
            retries -= 1
            sleep(retry_wait)
        else:
            break
        finally:
            urlcleanup()


def check_md5(file):
    print(f"Checking {file}", flush=True)
    md5_file = join(dirname(file), snakemake.params.md5_name)
    try:
        md5_df = pd.read_csv(
            md5_file,
            sep=r"\s+",
            header=None,
            names=["md5", "name"],
            engine="c",
            low_memory=False,
        )
        md5_df["name"] = md5_df["name"].apply(lambda n: basename(n))
        actual_md5 = md5_df["md5"][md5_df["name"] == basename(file)].values[0]
        file_md5 = md5(Path(file).read_bytes()).hexdigest()
    except Exception as e:
        remove(file)
        print(f"Error: {md5_file}: {e}", flush=True)
    else:
        if file_md5 != actual_md5:
            remove(file)
            print(f"Error: {file}: {file_md5} != {actual_md5}", flush=True)


print("\nDownloading NCBI genome assembly files")

summary_df = pd.read_csv(
    snakemake.input[0],
    sep="\t",
    index_col="assembly_accession",
    engine="c",
    low_memory=False,
)

genome_names = []
file_urls, files = [], []
md5_file_urls, md5_files = [], []
for acc in summary_df.index:
    ftp_dir_url = summary_df.loc[acc]["ftp_path"]
    genome_name = basename(urlparse(ftp_dir_url).path)
    genome_names.append(genome_name)
    genome_dir = join(snakemake.output[0], genome_name)
    for file_ext in snakemake.params.file_exts:
        file_name = "_".join([genome_name, file_ext])
        file_url = join(ftp_dir_url, file_name)
        file_urls.append(file_url)
        file = join(genome_dir, file_name)
        files.append(file)
    md5_file_url = join(ftp_dir_url, snakemake.params.md5_name)
    md5_file_urls.append(md5_file_url)
    md5_file = join(genome_dir, snakemake.params.md5_name)
    md5_files.append(md5_file)

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(
        url, file, snakemake.params.retries, snakemake.params.retry_wait
    )
    for url, file in zip(file_urls, files)
)

print("\nDownloading NCBI genome assembly md5sum files")

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(
        url, file, snakemake.params.retries, snakemake.params.retry_wait
    )
    for url, file in zip(md5_file_urls, md5_files)
)

print("\nChecking NCBI genome assembly md5sum files")

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(delayed(check_md5)(file) for file in [f for f in files if exists(f)])

# remove genome dirs with missing genomic sequence and annotation files
# testing if URLs exist (via request HEAD) before downloading is very slow so
# faster to download everything and then check/remove incomplete genomes
print("\nChecking for incomplete genomes")
for genome_name in genome_names:
    genome_dir = join(snakemake.output[0], genome_name)
    genome_file_names = listdir(genome_dir)
    if not any(
        n.endswith(tuple(snakemake.params.file_exts)) for n in genome_file_names
    ):
        print(f"Removing incomplete {genome_name}")
        rmtree(genome_dir, ignore_errors=True)
