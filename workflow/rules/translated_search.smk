rule unmapped_bam:
    input:
        PUFFERFISH_CDS_BAM_FILE,
    params:
        extra="--require-flags 12 --exclude-flags 256",
    output:
        PUFFERFISH_CDS_UNMAPPED_BAM_FILE,
    log:
        PUFFERFISH_CDS_UNMAPPED_BAM_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule unmapped_fastq:
    input:
        PUFFERFISH_CDS_UNMAPPED_BAM_FILE,
    output:
        PUFFERFISH_CDS_UNMAPPED_FASTQ_FILE,
    log:
        PUFFERFISH_CDS_UNMAPPED_FASTQ_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_FASTQ_WRAPPER
