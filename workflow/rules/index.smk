rule create_pufferfish_index:
    input:
        ref=NCBI_REFERENCE_FASTA_FILE,
    params:
        pufferfish=config["pufferfish"]["binary"],
        extra="--keepDuplicates",
    output:
        directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    resources:
        tmpdir=TEMP_DIR,
    threads: config["pufferfish"]["index"]["threads"]
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
