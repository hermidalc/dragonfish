from os.path import abspath


rule pufferfish_index:
    input:
        ref=REF_DEDUPED_ID_WITH_DECOY_FASTA_FILE,
        # decoys=GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE,
    params:
        pufferfish=abspath(join(config["pufferfish"]["bin_dir"], "pufferfish")),
        extra=config["pufferfish"]["index"]["extra"],
    output:
        directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    threads: PUFFERFISH_INDEX_THREADS
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
