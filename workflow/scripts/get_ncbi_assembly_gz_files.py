from os import makedirs
from os.path import join
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlcleanup, urlretrieve

import pandas as pd
from joblib import delayed, Parallel


def download_file(url, file, debug):
    if debug:
        print(url, flush=True)
    try:
        urlretrieve(url, filename=file)
        urlcleanup()
    except URLError as e:
        print(f"Missing {file}")


proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

summary_df = pd.read_csv(
    snakemake.input.summary_file,
    sep="\t",
    index_col="assembly_accession",
    keep_default_na=False,
)

file_urls, files = [], []
for acc in proteome_df["Genome assembly ID"].dropna():
    if acc in summary_df.index:
        dir_url = summary_df.loc[acc]["ftp_path"]
        dir_name = Path(dir_url).name
        for file_ext in snakemake.params.file_exts:
            file_name = "_".join([dir_name, file_ext])
            file_urls.append(join(dir_url, file_name))
            files.append(join(snakemake.output[0], file_name))
    else:
        print(f"Missing {acc}", flush=True)

makedirs(snakemake.output[0], mode=0o755, exist_ok=True)

Parallel(
    n_jobs=snakemake.threads,
    backend=snakemake.params.parallel_backend,
    verbose=snakemake.params.verbosity,
)(
    delayed(download_file)(url, file, snakemake.params.debug)
    for url, file in zip(file_urls, files)
)
