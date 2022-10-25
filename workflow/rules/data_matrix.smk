rule cds_bam_count_matrix:
    input:
        expand(PUFFERFISH_FILTERED_CDS_BAM_READ_QUANT_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
        data_col=1,
        collapse_techreps=True,
    output:
        PUFFERFISH_FILTERED_CDS_BAM_COUNT_MATRIX_FILE,
    log:
        PUFFERFISH_FILTERED_CDS_BAM_COUNT_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule cds_bam_count_eset:
    input:
        assay=PUFFERFISH_FILTERED_CDS_BAM_COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        PUFFERFISH_FILTERED_CDS_BAM_COUNT_ESET_FILE,
    log:
        PUFFERFISH_FILTERED_CDS_BAM_COUNT_ESET_LOG,
    wrapper:
        ESET_WRAPPER


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
