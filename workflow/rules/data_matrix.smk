rule cds_count_matrix:
    input:
        expand(PUFFERFISH_FILTERED_CDS_READ_QUANT_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
        data_col=1,
        collapse_techreps=True,
    output:
        PUFFERFISH_FILTERED_CDS_COUNT_MATRIX_FILE,
    log:
        PUFFERFISH_FILTERED_CDS_COUNT_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule uniprot_kb_dbxref_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        idmap=UNIPROT_KB_GENBANK_IDMAP_DBXREF_FILE,
        counts=PUFFERFISH_FILTERED_CDS_COUNT_MATRIX_FILE,
    output:
        UNIPROT_KB_DBXREF_COUNT_MATRIX_FILE,
    log:
        UNIPROT_KB_DBXREF_COUNT_MATRIX_LOG,
    threads: UNIPROT_KB_GENBANK_IDMAP_DBXREF_THREADS
    script:
        "../scripts/create_uniprot_kb_dbxref_count_matrix.py"


rule cedar_count_matrix:
    input:
        expand(CEDAR_READ_QUANT_FILE, zip, **EXPAND_PARAMS),
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
