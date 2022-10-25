from os.path import abspath


rule sam_to_bam:
    input:
        PUFFERFISH_MERGED_GENOMIC_CDS_FILTERED_SAM_FILE,
    params:
        extra="--bam",
    output:
        PUFFERFISH_MERGED_GENOMIC_CDS_FILTERED_BAM_FILE,
    log:
        PUFFERFISH_MERGED_GENOMIC_CDS_FILTERED_BAM_LOG,
    threads: SAMTOOLS_VIEW_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule filter_bam:
    input:
        bam=PUFFERFISH_MERGED_GENOMIC_CDS_FILTERED_BAM_FILE,
        qname=REF_CDS_FROM_GENOMIC_DEDUPED_ID_FILE,
    output:
        PUFFERFISH_CDS_FILTERED_BAM_FILE,
    log:
        PUFFERFISH_CDS_FILTERED_BAM_LOG,
    threads: SAMTOOLS_VIEW_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule bam_cds_read_quant:
    input:
        PUFFERFISH_CDS_FILTERED_BAM_FILE,
    params:
        read_count=True,
    output:
        PUFFERFISH_CDS_FILTERED_BAM_READ_QUANT_FILE,
    log:
        PUFFERFISH_CDS_FILTERED_BAM_READ_QUANT_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_BAM_WRAPPER


rule cedar_read_quant:
    input:
        PUFFERFISH_GENOMIC_PAM_FILE,
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
