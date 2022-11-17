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
        "../scripts/tsv_to_hdf.py"


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


rule cedar_count_matrix:
    conda:
        "../envs/pandas.yaml"
    input:
        CEDAR_READ_QUANT_FILE,
    output:
        CEDAR_COUNT_MATRIX_FILE,
    log:
        CEDAR_COUNT_MATRIX_LOG,
    script:
        "../scripts/cedar_count_matrix.py"


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
