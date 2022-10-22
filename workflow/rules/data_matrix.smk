rule cedar_count_matrix:
    input:
        expand(CEDAR_QUANT_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
        collapse_techreps=True,
    output:
        COUNT_MATRIX_FILE,
    log:
        COUNT_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule cedar_tpm_matrix:
    input:
        expand(CEDAR_QUANT_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
    output:
        TPM_MATRIX_FILE,
    log:
        TPM_MATRIX_LOG,
    wrapper:
        DATA_MATRIX_WRAPPER


rule count_eset:
    input:
        assay=COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        COUNT_ESET_FILE,
    log:
        COUNT_ESET_LOG,
    wrapper:
        ESET_WRAPPER


rule tpm_eset:
    input:
        assay=TPM_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        TPM_ESET_FILE,
    log:
        TPM_ESET_LOG,
    wrapper:
        ESET_WRAPPER
