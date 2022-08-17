rule create_pufferfish_index:
    input:
        ref=NCBI_REFERENCE_FASTA_FILE,
    params:
        pufferfish=config["pufferfish"]["binary"],
        extra="--keepDuplicates",
    output:
        tmp_dir=PUFFERFISH_TEMP_DIR,
        out_dir=directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    threads: config["pufferfish"]["index"]["threads"]
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
