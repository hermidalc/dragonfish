rule get_ncbi_acc2taxid_gz:
    params:
        NCBI_ACC2TAXID_GZ_URL,
    output:
        NCBI_ACC2TAXID_GZ_FILE,
    log:
        NCBI_ACC2TAXID_GZ_LOG,
    script:
        "../scripts/get_url_file.py"


rule gunzip_acc2taxid:
    input:
        NCBI_ACC2TAXID_GZ_FILE,
    output:
        NCBI_ACC2TAXID_FILE,
    log:
        NCBI_ACC2TAXID_LOG,
    script:
        "../scripts/unzip_file.py"


rule get_ncbi_taxdump_zip:
    params:
        NCBI_TAXDUMP_ZIP_URL,
    output:
        NCBI_TAXDUMP_ZIP_FILE,
    log:
        NCBI_TAXDUMP_ZIP_LOG,
    script:
        "../scripts/get_url_file.py"


rule unzip_ncbi_taxdump:
    input:
        NCBI_TAXDUMP_ZIP_FILE,
    output:
        NCBI_TAXDUMP_FILES,
    log:
        NCBI_TAXDUMP_FILES_LOG,
    script:
        "../scripts/unzip_file.py"
