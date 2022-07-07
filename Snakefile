import re
from glob import glob
from os import getcwd, mkdir, remove, walk
from os.path import exists, isdir, join, splitext
from shutil import rmtree

CONFIG_DIR = "config"
DATA_DIR = "data"
LOG_DIR = "logs"
RESULTS_DIR = "results"
RULES_DIR = "rules"
SCRIPTS_DIR = "scripts"


configfile: join(CONFIG_DIR, "config.yaml")


GENBANK_ASSEMBLY_DIR = join(DATA_DIR, "assembly")

UNIPROT_PROTEOME_METADATA_URL = config["uniprot_proteome_metadata_url"]
GENBANK_ASSEMBLY_REPORTS_URL = config["genbank_assembly_reports_url"]
GENBANK_ASSEMBLY_SUMMARY_FILENAMES = config["genbank_assembly_summary_filenames"]
GENBANK_ASSEMBLY_FILE_EXTS = config["genbank_assembly_file_exts"]

GENBANK_ASSEMBLY_SUMMARY_BASENAME = [
    splitext(filename)[0] for filename in GENBANK_ASSEMBLY_SUMMARY_FILENAMES
]

NCBI_API_EMAIL = config["ncbi_api_email"]
NCBI_API_KEY = config["ncbi_api_key"]

UNIPROT_PROTEOME_METADATA_FILE = join(DATA_DIR, "uniprot_proteome_metadata.tsv")
GENBANK_ASSEMBLY_SUMMARY_FILE = join(DATA_DIR, "{file_basename}.txt")
GENBANK_ASSEMBLY_SUMMARY_FILES = [
    join(DATA_DIR, filename) for filename in GENBANK_ASSEMBLY_SUMMARY_FILENAMES
]
GENBANK_MERGED_ASSEMBLY_SUMMARY_FILE = join(
    DATA_DIR, "assembly_summary_genbank_merged.txt"
)

UNIPROT_PROTEOME_METADATA_LOG = join(LOG_DIR, "get_uniprot_proteome_metadata.log")
GENBANK_ASSEMBLY_SUMMARY_LOG = join(LOG_DIR, "get_{file_basename}.log")
GENBANK_MERGED_ASSEMBLY_SUMMARY_LOG = join(
    LOG_DIR,
    "genbank_assembly_summaries.log",
)
GENBANK_ASSEMBLY_FILES_LOG = join(LOG_DIR, "get_genbank_assembly_files.log")


if not exists(LOG_DIR):
    mkdir(LOG_DIR, mode=0o755)


rule all:
    input:
        UNIPROT_PROTEOME_METADATA_FILE,
        expand(
            GENBANK_ASSEMBLY_SUMMARY_FILE,
            file_basename=GENBANK_ASSEMBLY_SUMMARY_BASENAME,
        ),
        GENBANK_MERGED_ASSEMBLY_SUMMARY_FILE,
        "data/finished",


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
        api_url=UNIPROT_PROTEOME_METADATA_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        UNIPROT_PROTEOME_METADATA_FILE,
    log:
        UNIPROT_PROTEOME_METADATA_LOG,
    shell:
        """
        python {params.scripts_dir}/get_uniprot_proteome_metadata.py \
        --api-url '{params.api_url}' \
        --out-file {output} \
        &> {log}
        """


rule get_genbank_assembly_summary:
    params:
        file_basename="{file_basename}",
        reports_url=GENBANK_ASSEMBLY_REPORTS_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        GENBANK_ASSEMBLY_SUMMARY_FILE,
    log:
        GENBANK_ASSEMBLY_SUMMARY_LOG,
    shell:
        """
        python {params.scripts_dir}/get_genbank_assembly_summary.py \
        --reports-url '{params.reports_url}' \
        --out-file {output} \
        &> {log}
        """


rule merge_genbank_assembly_summaries:
    input:
        GENBANK_ASSEMBLY_SUMMARY_FILES,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        GENBANK_MERGED_ASSEMBLY_SUMMARY_FILE,
    log:
        GENBANK_MERGED_ASSEMBLY_SUMMARY_LOG,
    shell:
        """
        python {params.scripts_dir}/merge_genbank_assembly_summaries.py \
        --summary-files {input} \
        --out-file {output} \
        &> {log}
        """


checkpoint get_genbank_assembly_files:
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=GENBANK_MERGED_ASSEMBLY_SUMMARY_FILE,
    params:
        file_exts=GENBANK_ASSEMBLY_FILE_EXTS,
        out_dir=GENBANK_ASSEMBLY_DIR,
        scripts_dir=SCRIPTS_DIR,
    output:
        directory(GENBANK_ASSEMBLY_DIR),
    log:
        GENBANK_ASSEMBLY_FILES_LOG,
    threads: 8
    shell:
        """
        python {params.scripts_dir}/get_genbank_assembly_files.py \
        --proteome-file {input.proteome_file} \
        --summary-file {input.summary_file} \
        --file-exts {params.file_exts} \
        --out-dir {params.out_dir} \
        --n-jobs {threads} \
        &> {log}
        """


def aggregate_genbank_assembly_files(wildcards):
    files = checkpoints.get_genbank_assembly_files.get(**wildcards).output[0]
    file_basenames, file_exts = glob_wildcards(join(files, "{i}.{ext}"))
    return expand(
        f"{GENBANK_ASSEMBLY_DIR}/{{i}}.{{ext}}", zip, file_basenames, file_exts
    )


rule finish:
    input:
        aggregate_genbank_assembly_files,
    output:
        touch("data/finished"),
    shell:
        "echo FINISHED"
