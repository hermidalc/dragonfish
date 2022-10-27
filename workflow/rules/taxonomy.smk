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
        "../scripts/get_url_file.py"


rule ncbi_acc2taxid_merged:
    input:
        expand(NCBI_ACC2TAXID_FILE, zip, **EXPAND_PARAMS),
    output:
        NCBI_ACC2TAXID_MERGED_FILE,
    log:
        NCBI_ACC2TAXID_MERGED_LOG,
    script:
        "../scripts/merge_ncbi_acc2taxid.py"


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
        "../scripts/get_url_file.py"


rule ncbi_taxdump_files:
    input:
        NCBI_TAXDUMP_ZIP_FILE,
    output:
        NCBI_TAXDUMP_FILES,
    log:
        NCBI_TAXDUMP_FILES_LOG,
    script:
        "../scripts/unzip_file.py"
