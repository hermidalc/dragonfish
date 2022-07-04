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
UNIPROT_PROTEOME_METADATA_FILE = join(DATA_DIR, "uniprot_proteome_metadata.pkl")
NCBI_FULL_GCA_NAME_FILE = join(DATA_DIR, "nci_full_gca_names.pkl")

UNIPROT_PROTEOME_METADATA_LOG = join(LOG_DIR, "get_uniprot_proteome_metadata.log")
NCBI_FULL_GCA_NAME_LOG = join(LOG_DIR, "get_ncbi_full_gca_names.log")


if not exists(LOG_DIR):
    mkdir(LOG_DIR, mode=0o755)


rule all:
    input:
        UNIPROT_PROTEOME_METADATA_FILE,
        NCBI_FULL_GCA_NAME_FILE,


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


rule get_ncbi_full_gca_names:
    input:
        UNIPROT_PROTEOME_METADATA_FILE,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_FULL_GCA_NAME_FILE,
    log:
        NCBI_FULL_GCA_NAME_LOG,
    shell:
        """
        python {params.scripts_dir}/get_ncbi_full_gca_names.py \
        --meta-file '{input}' \
        --name-file {output} \
        &> {log}
        """


# def input_for_get_ncbi_genome_file(wildcards):
#     metadata_df = load(checkpoints.get_uniprot_proteome_metadata.get().output)
#     conduit = Conduit("hermidalc@pitt.edu")
#     for acc in metadata_df["Genome assembly ID"]:
#         pipe = conduit.new_pipeline()
#         pid = pipe.add_search(
#             {"db": "assembly", "term": acc, "field": "Assembly Accession"}
#         )
#         pid = pipe.add_summary(dependency=pid)
#         analyzer = conduit.run(pipe)
#         results_size = analyzer.get_result().size()
#         if results_size > 1:
#             raise ValueError(f"{acc} query returned {results_size:d} results")
#         ftp_dir_url = list(analyzer.get_result().summaries.values())[0][
#             "ftppath_genbank"
#         ]
#         url_parts = urlparse(ftp_dir_url)
#         ftp = FTP(url_parts.hostname)
#         ftp.login()
#         gca_name = split(url_parts.path)[-1]
#         get_filenames = [
#             f"{gca_name}_cds_from_genomic.fna.gz",
#             f"{gca_name}_genomic.fna.gz",
#             f"{gca_name}_genomic.gtf.gz",
#             "md5checksums.txt",
#         ]
#         ftp_dir_filenames = ftp.nlst(url_parts.path)
#         for filename in get_filenames:
#             if filename not in ftp_dir_filenames:
#                 raise ValueError(f"{filename} not in {ftp_dir_url}")
#         ftp_files = [f"{ftp_dir_url}/{filename}" for filename in get_filenames]
#         break
#     return ftp_files


# rule get_ncbi_genome_file:
#     input:
#         input_for_get_ncbi_genome_file,
#     output:
#         UNIPROT_PROTEOME_METADATA_FILE,
#     log:
#         UNIPROT_PROTEOME_METADATA_LOG,
#     shell:
#         """
#         python {params.scripts_dir}/get_ncbi_genome_file.py \
#         --ftp-file-url '{params.api_url}' \
#         --out-file {output} \
#         &> {log}
#         """
