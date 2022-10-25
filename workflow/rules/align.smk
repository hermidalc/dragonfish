from os.path import abspath


rule pufferfish_align_pe:
    input:
        unpack(lambda wc: get_fq(wc, trimmed=True)),
        index=PUFFERFISH_INDEX_DIR,
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
