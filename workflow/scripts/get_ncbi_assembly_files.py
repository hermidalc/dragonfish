from hashlib import md5
from os import listdir
from os.path import basename, dirname, exists, join
from pathlib import Path
from shutil import rmtree
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import urlcleanup, urlretrieve

import numpy as np
import pandas as pd
from joblib import delayed, Parallel
from snakemake.utils import makedirs

md5_file_name = "md5checksums.txt"


def download_file(url, file, debug):
    if debug:
        print(f"Downloading {url}", flush=True)
    makedirs(dirname(file))
    try:
        urlretrieve(url, filename=file)
    except URLError as e:
        print(f"Skipped: {url}: {e}", flush=True)
    else:
        urlcleanup()


def check_md5(file, debug):
    md5_file = join(dirname(file), md5_file_name)
    if debug:
        print(f"Checking {md5_file}", flush=True)
    try:
        md5_df = pd.read_csv(
            md5_file, delim_whitespace=True, header=None, names=["md5", "name"]
        )
        md5_df["name"] = md5_df["name"].apply(lambda n: basename(n))
        actual_md5 = md5_df["md5"][md5_df["name"] == basename(file)].values[0]
        file_md5 = md5(Path(file).read_bytes()).hexdigest()
        if file_md5 != actual_md5:
            print(f"{basename(file)} md5 {file_md5} != {actual_md5}")
    except Exception as e:
        print(f"MD5 Error: {md5_file}: {e}")


print(
    "\nGetting NCBI genomes to download from UniProt Proteomes and NCBI assembly metadata"
)

summary_df = pd.read_csv(
    snakemake.input.summary_file,
    sep="\t",
    index_col="assembly_accession",
    low_memory=False,
)
proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

genome_names = []
file_urls, files = [], []
md5_file_urls, md5_files = [], []
for acc in proteome_df["Genome assembly ID"]:
    if acc in summary_df.index:
        ftp_dir_url = summary_df.loc[acc]["ftp_path"]
        if pd.notna(ftp_dir_url):
            genome_name = basename(urlparse(ftp_dir_url).path)
            genome_names.append(genome_name)
            # genome_dir_name = genome_name.removesuffix(".gz").removesuffix(".GZ")
            genome_dir_name = genome_name
            genome_dir = join(snakemake.output[0], genome_dir_name)
            for file_ext in snakemake.params.file_exts:
                file_url = join(ftp_dir_url, "_".join([genome_name, file_ext]))
                file_urls.append(file_url)
                file_name = "_".join([genome_dir_name, file_ext])
                file = join(genome_dir, file_name)
                files.append(file)
            md5_file_url = join(ftp_dir_url, md5_file_name)
            md5_file_urls.append(md5_file_url)
            md5_file = join(genome_dir, md5_file_name)
            md5_files.append(md5_file)
        else:
            print(f"No URL {acc}")
    else:
        print(f"Missing {acc}")

genome_name_series = pd.Series(genome_names)
dup_genome_names = genome_name_series[genome_name_series.duplicated()].values
assert (
    dup_genome_names.size == 0
), "Duplicate genome names found in UniProt Proteomes:\n{}".format(
    "\n".join(dup_genome_names)
)

print("\nDownloading NCBI genome assembly files")

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(url, file, snakemake.params.debug)
    for url, file in zip(file_urls, files)
)

print("\nDownloading NCBI genome assembly md5sum files")

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(url, file, snakemake.params.debug)
    for url, file in zip(md5_file_urls, md5_files)
)

print("\nChecking NCBI genome assembly md5sum files")

downloaded_files = [f for f in files if exists(f)]

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.backend,
    verbose=snakemake.params.verbosity,
)(delayed(check_md5)(file, snakemake.params.debug) for file in downloaded_files)

# XXX: not sure yet if to include genomes with missing CDSs or GTF/GFF
# remove incomplete file groups
# testing if URLs exist (via request HEAD) before downloading is very slow so actually
# faster to download everything and then check/remove incomplete file groups
print("\nChecking for incomplete genome files")
for genome_name in genome_names:
    # genome_dir_name = genome_name.removesuffix(".gz").removesuffix(".GZ")
    genome_dir_name = genome_name
    genome_dir = join(snakemake.output[0], genome_dir_name)
    invalid_file_group = False
    genome_file_names = listdir(genome_dir)
    for file_ext in snakemake.params.file_exts:
        if not np.any([n.endswith(file_ext) for n in genome_file_names]):
            invalid_file_group = True
            break
    if invalid_file_group:
        print(f"Incomplete {genome_dir}")
        rmtree(genome_dir)
