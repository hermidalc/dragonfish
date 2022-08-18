rule create_pufferfish_index:
    input:
        ref=NCBI_REFERENCE_FASTA_FILE,
    params:
        pufferfish=config["pufferfish"]["binary"],
        extra="--keepDuplicates",
        tmp_dir=TEMP_DIR,
    output:
        out_dir=directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    threads: config["pufferfish"]["index"]["threads"]
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
