import re
from ftplib import FTP
from glob import glob
from os import getcwd, mkdir, remove, walk
from os.path import exists, isdir, join, split
from shutil import rmtree
from urllib.parse import urlparse

import pandas as pd
from joblib import load
from entrezpy.conduit import Conduit


CONFIG_DIR = "config"
DATA_DIR = "data"
LOG_DIR = "logs"
RESULTS_DIR = "results"
RULES_DIR = "rules"
SCRIPTS_DIR = "scripts"


configfile: join(CONFIG_DIR, "config.yaml")


UNIPROT_PROTEOME_METADATA_URL = config["uniprot_proteome_metadata_url"]
UNIPROT_PROTEOME_METADATA_FILE = join(DATA_DIR, "uniprot_proteome_metadata.tsv")
NCBI_GCA_NAME_FILE = join(DATA_DIR, "ncbi_gca_names.tsv")

UNIPROT_PROTEOME_METADATA_LOG = join(LOG_DIR, "get_uniprot_proteome_metadata.log")
NCBI_GCA_NAME_LOG = join(LOG_DIR, "get_ncbi_gca_names.log")


if not exists(LOG_DIR):
    mkdir(LOG_DIR, mode=0o755)


rule all:
    input:
        UNIPROT_PROTEOME_METADATA_FILE,
        NCBI_GCA_NAME_FILE,


rule clean:
    run:
        for file in glob(join(DATA_DIR, "*")):
            if isdir(file):
                rmtree(file)
            else:
                remove(file)
        for file in glob(join(LOG_DIR, "*")):
            remove(file)
        for dirpath, dirnames, filenames in sorted(walk(getcwd())):
            for dirname in dirnames:
                if dirname == "__pycache__":
                    rmtree(join(dirpath, dirname))


rule get_uniprot_proteome_metadata:
    params:
        scripts_dir=SCRIPTS_DIR,
        api_url=UNIPROT_PROTEOME_METADATA_URL,
    output:
        UNIPROT_PROTEOME_METADATA_FILE,
    log:
        UNIPROT_PROTEOME_METADATA_LOG,
    shell:
        """
        python {params.scripts_dir}/get_uniprot_proteome_metadata.py \
        --api-url '{params.api_url}' \
        --meta-file {output} \
        &> {log}
        """


rule get_ncbi_gca_names:
    input:
        UNIPROT_PROTEOME_METADATA_FILE,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_GCA_NAME_FILE,
    log:
        NCBI_GCA_NAME_LOG,
    threads: 5
    shell:
        """
        python {params.scripts_dir}/get_ncbi_gca_names.py \
        --meta-file {input} \
        --name-file {output} \
        --n-jobs {threads}
        &> {log}
        """
