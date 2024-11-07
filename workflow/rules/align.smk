from os.path import abspath


rule pufferfish_align:
    input:
        unpack(lambda wc: get_fq(wc, data_type="filtered")),
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


rule pufferfish_unmapped_bam:
    input:
        PUFFERFISH_SAM_FILE,
    params:
        extra="--require-flags 12 --exclude-flags 256",
    output:
        PUFFERFISH_UNMAPPED_BAM_FILE,
    log:
        PUFFERFISH_UNMAPPED_BAM_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule pufferfish_unmapped_fastq_pe:
    input:
        PUFFERFISH_UNMAPPED_BAM_FILE,
    output:
        PUFFERFISH_UNMAPPED_FASTQ_R1_FILE,
        PUFFERFISH_UNMAPPED_FASTQ_R2_FILE,
    log:
        PUFFERFISH_UNMAPPED_FASTQ_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_FASTQ_WRAPPER
