from os.path import abspath


rule cedar_read_quant:
    input:
        PUFFERFISH_ALIGN_FILE,
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
