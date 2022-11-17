from os.path import abspath


rule featurecounts_cds_read_quant:
    input:
        expand(PUFFERFISH_SAM_FILE, zip, **EXPAND_PARAMS),
    output:
        FEATURECOUNTS_CDS_READ_QUANT_FILE,
    log:
        FEATURECOUNTS_CDS_READ_QUANT_LOG,
    threads: FEATURECOUNTS_THREADS
    wrapper:
        FEATURECOUNTS_WRAPPER


rule cedar_read_quant:
    input:
        pam=PUFFERFISH_PAM_FILE,
        taxtree=NCBI_TAXDUMP_FIXED_NODE_FILE,
        ref2tax=REF_ID2TAXID_FILE,
    params:
        cedar=abspath(join(config["pufferfish"]["bin_dir"], "cedar")),
        extra=config["pufferfish"]["cedar"]["extra"],
    output:
        CEDAR_READ_QUANT_FILE,
    log:
        CEDAR_READ_QUANT_LOG,
    threads: CEDAR_READ_QUANT_THREADS
    wrapper:
        CEDAR_READ_QUANT_WRAPPER
