rule featurecounts_cds_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        FEATURECOUNTS_CDS_READ_QUANT_FILE,
    output:
        FEATURECOUNTS_CDS_COUNT_MATRIX_FILE,
    log:
        FEATURECOUNTS_CDS_COUNT_MATRIX_LOG,
    script:
        "../scripts/featurecounts_count_matrix.py"


rule featurecounts_dbxref_cds_count_matrix:
    conda:
        "../envs/vaex.yaml"
    input:
        counts=FEATURECOUNTS_CDS_COUNT_MATRIX_FILE,
        idmap=UNIPROT_KB_GENBANK_IDMAP_DBXREF_FILE,
    output:
        FEATURECOUNTS_DBXREF_CDS_COUNT_MATRIX_FILE,
    log:
        FEATURECOUNTS_DBXREF_CDS_COUNT_MATRIX_LOG,
    threads: UNIPROT_KB_GENBANK_IDMAP_DBXREF_THREADS
    script:
        "../scripts/uniprot_kb_dbxref_count_matrix.py"


rule cedar_tax_count_matrix:
    conda:
        "../envs/pandas.yaml"
    input:
        expand(CEDAR_TAX_READ_QUANT_FILE, zip, **EXPAND_PARAMS),
    output:
        CEDAR_TAX_COUNT_MATRIX_FILE,
    log:
        CEDAR_TAX_COUNT_MATRIX_LOG,
    script:
        "../scripts/cedar_tax_count_matrix.py"


rule cedar_tax_count_eset:
    input:
        assay=CEDAR_TAX_COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        CEDAR_TAX_COUNT_ESET_FILE,
    log:
        CEDAR_TAX_COUNT_ESET_LOG,
    wrapper:
        ESET_WRAPPER
