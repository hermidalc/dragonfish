rule ncbi_acc2taxid:
    params:
        NCBI_ACC2TAXID_URL,
    output:
        NCBI_ACC2TAXID_FILE,
    log:
        NCBI_ACC2TAXID_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"


rule ncbi_acc2taxid_filtered_merged:
    conda:
        "../envs/pandas.yaml"
    input:
        files=expand(NCBI_ACC2TAXID_FILE, zip, **EXPAND_PARAMS),
        proteomes=UNIPROT_PROTEOMES_FILE,
    output:
        NCBI_ACC2TAXID_FILTERED_MERGED_FILE,
    log:
        NCBI_ACC2TAXID_FILTERED_MERGED_LOG,
    script:
        "../scripts/ncbi_acc2taxid_filtered_merged.py"


rule ncbi_taxdump_zip:
    params:
        NCBI_TAXDUMP_ZIP_URL,
    output:
        NCBI_TAXDUMP_ZIP_FILE,
    log:
        NCBI_TAXDUMP_ZIP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"


rule ncbi_taxdump_files:
    input:
        NCBI_TAXDUMP_ZIP_FILE,
    output:
        NCBI_TAXDUMP_FILES,
    log:
        NCBI_TAXDUMP_FILES_LOG,
    script:
        "../scripts/unzip_file.py"
