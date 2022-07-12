import re
from glob import glob
from os import getcwd, mkdir, remove, walk
from os.path import basename, dirname, exists, isdir, join, splitext
from shutil import rmtree

CONFIG_DIR = "config"
DATA_DIR = "data"
LOG_DIR = "logs"
RESULTS_DIR = "results"
RULES_DIR = "rules"
SCRIPTS_DIR = "scripts"

NCBI_TAXONOMY_DIR = join(DATA_DIR, "taxonomy")
NCBI_GENOME_DIR = join(DATA_DIR, "genomes")
NCBI_ASSEMBLY_DIR = join(NCBI_GENOME_DIR, "assemblies")
UNIPROT_DIR = join(DATA_DIR, "uniprot")


configfile: join(CONFIG_DIR, "config.yaml")


NCBI_API_EMAIL = config["ncbi_api_email"]
NCBI_API_KEY = config["ncbi_api_key"]

NCBI_TAXONOMY_URL = config["ncbi_taxonomy_url"]
NCBI_TAXONOMY_FILE_PATHS = config["ncbi_taxonomy_file_paths"]
NCBI_ASSEMBLY_SUMMARY_URL = config["ncbi_assembly_summary_url"]
NCBI_ASSEMBLY_SUMMARY_FILENAMES = config["ncbi_assembly_summary_filenames"]
NCBI_ASSEMBLY_FILE_EXTS = config["ncbi_assembly_file_exts"]
UNIPROT_PROTEOME_METADATA_URL = config["uniprot_proteome_metadata_url"]

NCBI_TAXONOMY_FILE_DIRNAME, NCBI_TAXONOMY_FILE_BASENAME, NCBI_TAXONOMY_FILE_EXT = zip(
    *(
        (
            dirname(file_path),
            splitext(basename(file_path))[0],
            splitext(basename(file_path))[1].replace(".", "", 1),
        )
        for file_path in NCBI_TAXONOMY_FILE_PATHS
    )
)
NCBI_ASSEMBLY_SUMMARY_FILE_BASENAME, NCBI_ASSEMBLY_SUMMARY_FILE_EXT = zip(
    *(
        (splitext(filename)[0], splitext(filename)[1].replace(".", "", 1))
        for filename in NCBI_ASSEMBLY_SUMMARY_FILENAMES
    )
)

NCBI_TAXONOMY_FILE = join(
    NCBI_TAXONOMY_DIR, "{tax_dirname}", "{tax_basename}.{tax_ext}"
)
NCBI_TAXONOMY_FILE_URL = join(
    NCBI_TAXONOMY_URL, "{tax_dirname}", "{tax_basename}.{tax_ext}"
)
NCBI_ASSEMBLY_SUMMARY_FILE = join(NCBI_GENOME_DIR, "{asm_basename}.{asm_ext}")
NCBI_ASSEMBLY_SUMMARY_FILE_URL = join(
    NCBI_ASSEMBLY_SUMMARY_URL, "{asm_basename}.{asm_ext}"
)
NCBI_ASSEMBLY_SUMMARY_FILES = [
    join(NCBI_GENOME_DIR, filename) for filename in NCBI_ASSEMBLY_SUMMARY_FILENAMES
]
NCBI_MERGED_ASSEMBLY_SUMMARY_FILE = join(NCBI_GENOME_DIR, "assembly_summary_merged.txt")
UNIPROT_PROTEOME_METADATA_FILE = join(UNIPROT_DIR, "uniprot_proteome_metadata.tsv")

NCBI_TAXONOMY_LOG = join(LOG_DIR, "get_{tax_dirname}_{tax_basename}_{tax_ext}.log")
NCBI_ASSEMBLY_SUMMARY_LOG = join(LOG_DIR, "get_{asm_basename}_{asm_ext}.log")
NCBI_MERGED_ASSEMBLY_SUMMARY_LOG = join(LOG_DIR, "merge_ncbi_assembly_summaries.log")
UNIPROT_PROTEOME_METADATA_LOG = join(LOG_DIR, "get_uniprot_proteome_metadata.log")
NCBI_ASSEMBLY_FILE_LOG = join(LOG_DIR, "get_ncbi_assembly_files.log")


if not exists(LOG_DIR):
    mkdir(LOG_DIR, mode=0o755)


rule all:
    input:
        expand(
            NCBI_TAXONOMY_FILE,
            zip,
            tax_dirname=NCBI_TAXONOMY_FILE_DIRNAME,
            tax_basename=NCBI_TAXONOMY_FILE_BASENAME,
            tax_ext=NCBI_TAXONOMY_FILE_EXT,
        ),
        expand(
            NCBI_ASSEMBLY_SUMMARY_FILE,
            zip,
            asm_basename=NCBI_ASSEMBLY_SUMMARY_FILE_BASENAME,
            asm_ext=NCBI_ASSEMBLY_SUMMARY_FILE_EXT,
        ),
        UNIPROT_PROTEOME_METADATA_FILE,
        NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,


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


rule get_ncbi_taxonomy:
    params:
        file_url=NCBI_TAXONOMY_FILE_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_TAXONOMY_FILE,
    log:
        NCBI_TAXONOMY_LOG,
    shell:
        """
        python {params.scripts_dir}/get_url_file.py \
        --file-url '{params.file_url}' \
        --out-file {output} \
        &> {log}
        """


rule get_ncbi_assembly_summary:
    params:
        file_url=NCBI_ASSEMBLY_SUMMARY_FILE_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_ASSEMBLY_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_SUMMARY_LOG,
    shell:
        """
        python {params.scripts_dir}/get_url_file.py \
        --file-url '{params.file_url}' \
        --out-file {output} \
        &> {log}
        """


rule merge_ncbi_assembly_summaries:
    input:
        NCBI_ASSEMBLY_SUMMARY_FILES,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
    log:
        NCBI_MERGED_ASSEMBLY_SUMMARY_LOG,
    shell:
        """
        python {params.scripts_dir}/merge_ncbi_assembly_summaries.py \
        --summary-files {input} \
        --out-file {output} \
        &> {log}
        """


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


checkpoint get_ncbi_assembly_files:
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
    params:
        file_exts=NCBI_ASSEMBLY_FILE_EXTS,
        out_dir=NCBI_ASSEMBLY_DIR,
        scripts_dir=SCRIPTS_DIR,
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_FILE_LOG,
    threads: 8
    shell:
        """
        python {params.scripts_dir}/get_ncbi_assembly_files.py \
        --proteome-file {input.proteome_file} \
        --summary-file {input.summary_file} \
        --file-exts {params.file_exts} \
        --out-dir {params.out_dir} \
        --n-jobs {threads} \
        &> {log}
        """


def aggregate_ncbi_assembly_files(wildcards):
    files = checkpoints.get_ncbi_assembly_files.get(**wildcards).output[0]
    basenames, exts = glob_wildcards(join(files, "{asf_basename}.{asf_ext}"))
    return expand(
        f"{NCBI_ASSEMBLY_DIR}/{{asf_basename}}.{{asf_ext}}", zip, basenames, exts
    )


rule finish:
    input:
        aggregate_ncbi_assembly_files,
    output:
        touch("data/finished"),
    shell:
        "echo FINISHED"
