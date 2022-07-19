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

NCBI_TAXONOMY_DIR = join(DATA_DIR, "taxonomy")
NCBI_GENOME_DIR = join(DATA_DIR, "genomes")
NCBI_ASSEMBLY_DIR = join(NCBI_GENOME_DIR, "assemblies")
UNIPROT_DIR = join(DATA_DIR, "uniprot")


configfile: join(CONFIG_DIR, "config.yaml")


NCBI_API_EMAIL = config["ncbi"]["api_email"]
NCBI_API_KEY = config["ncbi"]["api_key"]

NCBI_ACC2TAXID_URL = config["ncbi"]["acc2taxid_url"]
NCBI_ACC2TAXID_GZ_FILENAMES = config["ncbi"]["acc2taxid_gz_filenames"]
NCBI_TAXDUMP_URL = config["ncbi"]["taxdump_url"]
NCBI_TAXDUMP_ZIP_FILENAME = config["ncbi"]["taxdump_zip_filename"]
NCBI_TAXDUMP_FILENAMES = config["ncbi"]["taxdump_filenames"]

NCBI_ASSEMBLY_SUMMARY_URL = config["ncbi"]["assembly_summary_url"]
NCBI_ASSEMBLY_SUMMARY_FILENAMES = config["ncbi"]["assembly_summary_filenames"]
NCBI_ASSEMBLY_GZ_FILE_EXTS = config["ncbi"]["assembly_gz_file_exts"]
UNIPROT_PROTEOME_METADATA_URL = config["uniprot"]["proteome_metadata_url"]

NCBI_ACC2TAXID_FILENAMES = [
    splitext(filename)[0] for filename in NCBI_ACC2TAXID_GZ_FILENAMES
]
NCBI_ACC2TAXID_GZ_URL = join(NCBI_ACC2TAXID_URL, "{a2t_filename}.gz")
NCBI_ACC2TAXID_GZ_FILE = join(NCBI_TAXONOMY_DIR, "{a2t_filename}.gz")
NCBI_ACC2TAXID_FILE = join(NCBI_TAXONOMY_DIR, "{a2t_filename}")
NCBI_TAXDUMP_ZIP_URL = join(NCBI_TAXDUMP_URL, NCBI_TAXDUMP_ZIP_FILENAME)
NCBI_TAXDUMP_ZIP_FILE = join(NCBI_TAXONOMY_DIR, NCBI_TAXDUMP_ZIP_FILENAME)
NCBI_TAXDUMP_FILES = [
    join(NCBI_TAXONOMY_DIR, filename) for filename in NCBI_TAXDUMP_FILENAMES
]
NCBI_ASSEMBLY_SUMMARY_BASENAMES, NCBI_ASSEMBLY_SUMMARY_EXTS = zip(
    *(
        (splitext(filename)[0], splitext(filename)[1].replace(".", "", 1))
        for filename in NCBI_ASSEMBLY_SUMMARY_FILENAMES
    )
)
NCBI_ASSEMBLY_SUMMARY_FILE_URL = join(
    NCBI_ASSEMBLY_SUMMARY_URL, "{asu_basename}.{asu_ext}"
)
NCBI_ASSEMBLY_SUMMARY_FILE = join(NCBI_GENOME_DIR, "{asu_basename}.{asu_ext}")
NCBI_ASSEMBLY_SUMMARY_FILES = [
    join(NCBI_GENOME_DIR, filename) for filename in NCBI_ASSEMBLY_SUMMARY_FILENAMES
]
NCBI_MERGED_ASSEMBLY_SUMMARY_FILE = join(NCBI_GENOME_DIR, "assembly_summary_merged.txt")
UNIPROT_PROTEOME_METADATA_FILE = join(UNIPROT_DIR, "uniprot_proteome_metadata.tsv")
NCBI_ASSEMBLY_GZ_FILE = join(NCBI_ASSEMBLY_DIR, "{asm_basename}.{asm_ext}.gz")
NCBI_ASSEMBLY_FILE = join(NCBI_ASSEMBLY_DIR, "{asm_basename}.{asm_ext}")
NCBI_ASSEMBLY_EXTS = [
    splitext(ext)[1].replace(".", "", 1) for ext in NCBI_ASSEMBLY_GZ_FILE_EXTS
]

NCBI_ACC2TAXID_GZ_LOG = join(LOG_DIR, "get_{a2t_filename}_gz.log")
NCBI_ACC2TAXID_LOG = join(LOG_DIR, "gunzip_{a2t_filename}_gz.log")
NCBI_TAXDUMP_ZIP_LOG = join(LOG_DIR, "get_ncbi_taxdump_zip.log")
NCBI_TAXDUMP_FILES_LOG = join(LOG_DIR, "unzip_ncbi_taxdump.log")
NCBI_ASSEMBLY_SUMMARY_LOG = join(LOG_DIR, "get_{asu_basename}_{asu_ext}.log")
NCBI_MERGED_ASSEMBLY_SUMMARY_LOG = join(LOG_DIR, "merge_ncbi_assembly_summaries.log")
UNIPROT_PROTEOME_METADATA_LOG = join(LOG_DIR, "get_uniprot_proteome_metadata.log")
NCBI_ASSEMBLY_GZ_FILES_LOG = join(LOG_DIR, "get_ncbi_assembly_gz_files.log")
NCBI_ASSEMBLY_FILE_LOG = join(LOG_DIR, "gunzip_{asm_basename}_{asm_ext}_gz.log")

if not exists(LOG_DIR):
    mkdir(LOG_DIR, mode=0o755)


wildcard_constraints:
    asu_basename="|".join(set(re.escape(NCBI_ASSEMBLY_SUMMARY_BASENAMES))),
    asu_ext="|".join(set(re.escape(NCBI_ASSEMBLY_SUMMARY_EXTS))),
    a2t_filename="|".join(set(re.escape(NCBI_ACC2TAXID_FILENAMES))),


rule all:
    input:
        expand(NCBI_ACC2TAXID_GZ_FILE, zip, a2t_filename=NCBI_ACC2TAXID_FILENAMES),
        expand(NCBI_ACC2TAXID_FILE, zip, a2t_filename=NCBI_ACC2TAXID_FILENAMES),
        NCBI_TAXDUMP_ZIP_FILE,
        NCBI_TAXDUMP_FILES,
        expand(
            NCBI_ASSEMBLY_SUMMARY_FILE,
            zip,
            asu_basename=NCBI_ASSEMBLY_SUMMARY_BASENAMES,
            asu_ext=NCBI_ASSEMBLY_SUMMARY_EXTS,
        ),
        NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
        UNIPROT_PROTEOME_METADATA_FILE,
        join(NCBI_GENOME_DIR, "aggregated.txt"),


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


rule get_ncbi_acc2taxid_gz:
    params:
        file_url=NCBI_ACC2TAXID_GZ_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_ACC2TAXID_GZ_FILE,
    log:
        NCBI_ACC2TAXID_GZ_LOG,
    shell:
        """
        python {params.scripts_dir}/get_url_file.py \
        --file-url '{params.file_url}' \
        --out-file {output} \
        1> {log} 2>&1
        """


rule gunzip_acc2taxid:
    input:
        NCBI_ACC2TAXID_GZ_FILE,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_ACC2TAXID_FILE,
    log:
        NCBI_ACC2TAXID_LOG,
    shell:
        """
        python {params.scripts_dir}/gunzip_file.py \
        --file {input} \
        1> {log} 2>&1
        """


rule get_ncbi_taxdump_zip:
    params:
        file_url=NCBI_TAXDUMP_ZIP_URL,
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_TAXDUMP_ZIP_FILE,
    log:
        NCBI_TAXDUMP_ZIP_LOG,
    shell:
        """
        python {params.scripts_dir}/get_url_file.py \
        --file-url '{params.file_url}' \
        --out-file {output} \
        1> {log} 2>&1
        """


rule unzip_ncbi_taxdump:
    input:
        NCBI_TAXDUMP_ZIP_FILE,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_TAXDUMP_FILES,
    log:
        NCBI_TAXDUMP_FILES_LOG,
    shell:
        """
        python {params.scripts_dir}/unzip_file.py \
        --file {input} \
        --members {output} \
        1> {log} 2>&1
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
        1> {log} 2>&1
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
        1> {log} 2>&1
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
        1> {log} 2>&1
        """


checkpoint get_ncbi_assembly_gz_files:
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
    params:
        file_exts=NCBI_ASSEMBLY_GZ_FILE_EXTS,
        out_dir=NCBI_ASSEMBLY_DIR,
        scripts_dir=SCRIPTS_DIR,
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_GZ_FILES_LOG,
    threads: 10
    shell:
        """
        python {params.scripts_dir}/get_ncbi_assembly_gz_files.py \
        --proteome-file {input.proteome_file} \
        --summary-file {input.summary_file} \
        --file-exts {params.file_exts} \
        --out-dir {params.out_dir} \
        --n-jobs {threads} \
        1> {log} 2>&1
        """


rule gunzip_ncbi_assembly_gz_file:
    input:
        NCBI_ASSEMBLY_GZ_FILE,
    params:
        scripts_dir=SCRIPTS_DIR,
    output:
        NCBI_ASSEMBLY_FILE,
    log:
        NCBI_ASSEMBLY_FILE_LOG,
    shell:
        """
        python {params.scripts_dir}/gunzip_file.py \
        --file {input} \
        1> {log} 2>&1
        """


def gather_ncbi_assembly_files(wildcards):
    out_dir = checkpoints.get_ncbi_assembly_gz_files.get(**wildcards).output[0]
    basenames, exts = glob_wildcards(join(out_dir, "{asm_basename}.{asm_ext}.gz"))
    return expand(
        f"{out_dir}/{{asm_basename}}.{{asm_ext}}",
        zip,
        asm_basename=basenames,
        asm_ext=exts,
    )


rule finished:
    input:
        gather_ncbi_assembly_files,
    output:
        join(NCBI_GENOME_DIR, "aggregated.txt"),
    shell:
        """
        echo {input} > {output}
        """
