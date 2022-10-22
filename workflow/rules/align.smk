rule pufferfish_align_pe:
    input:
        index=PUFFERFISH_INDEX_DIR,
        fq1=TRIMMED_FASTQ1_FILE,
        fq2=TRIMMED_FASTQ2_FILE,
    params:
        pufferfish=abspath(join(config["pufferfish"]["bin_dir"], "pufferfish")),
        extra=config["pufferfish"]["align"]["extra"],
    output:
        PUFFERFISH_ALIGN_FILE,
    log:
        PUFFERFISH_ALIGN_LOG,
    threads: PUFFERFISH_ALIGN_THREADS
    wrapper:
        PUFFERFISH_ALIGN_WRAPPER
