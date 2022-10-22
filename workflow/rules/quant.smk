rule quant_reads_align_cedar:
    input:
        PUFFERFISH_ALIGN_FILE,
    params:
        cedar=abspath(join(config["pufferfish"]["bin_dir"], "cedar")),
        extra=config["pufferfish"]["cedar"]["extra"],
    output:
        CEDAR_QUANT_FILE,
    log:
        CEDAR_QUANT_LOG,
    threads: CEDAR_QUANT_THREADS
    wrapper:
        CEDAR_QUANT_WRAPPER
