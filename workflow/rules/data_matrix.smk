rule cds_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        expand(PUFFERFISH_CDS_READ_QUANT_HDF_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
        data_col=1,
        collapse_techreps=True,
    output:
        PUFFERFISH_CDS_COUNT_MATRIX_FILE,
    log:
        PUFFERFISH_CDS_COUNT_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule cds_dbxref_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        counts=PUFFERFISH_CDS_COUNT_MATRIX_FILE,
        idmap=UNIPROT_KB_GENBANK_IDMAP_DBXREF_FILE,
    output:
        PUFFERFISH_CDS_DBXREF_COUNT_MATRIX_FILE,
    log:
        PUFFERFISH_CDS_DBXREF_COUNT_MATRIX_LOG,
    threads: UNIPROT_KB_GENBANK_IDMAP_DBXREF_THREADS
    script:
        "../scripts/uniprot_kb_dbxref_count_matrix.py"


rule cds_dbxref_count_eset:
    input:
        assay=PUFFERFISH_CDS_DBXREF_COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        PUFFERFISH_CDS_DBXREF_COUNT_ESET_FILE,
    log:
        PUFFERFISH_CDS_DBXREF_COUNT_ESET_LOG,
    wrapper:
        ESET_WRAPPER


rule cedar_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        expand(CEDAR_READ_QUANT_HDF_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
        data_col=4,
        collapse_techreps=True,
    output:
        CEDAR_COUNT_MATRIX_FILE,
    log:
        CEDAR_COUNT_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule cedar_count_eset:
    input:
        assay=CEDAR_COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        CEDAR_COUNT_ESET_FILE,
    log:
        CEDAR_COUNT_ESET_LOG,
    wrapper:
        ESET_WRAPPER
