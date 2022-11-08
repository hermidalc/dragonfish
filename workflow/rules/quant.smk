from os.path import abspath


rule cds_and_genomic_bam:
    input:
        PUFFERFISH_GENOMIC_AND_CDS_SAM_FILE,
    params:
        extra="--bam",
    output:
        PUFFERFISH_GENOMIC_AND_CDS_BAM_FILE,
    log:
        PUFFERFISH_GENOMIC_AND_CDS_BAM_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule cds_bam:
    input:
        bam=PUFFERFISH_GENOMIC_AND_CDS_BAM_FILE,
        qname=REF_ASSEMBLY_CDS_FROM_GENOMIC_DEDUPED_QNAME_FILE,
    output:
        PUFFERFISH_CDS_BAM_FILE,
    log:
        PUFFERFISH_CDS_BAM_LOG,
    threads: SAMTOOLS_THREADS
    wrapper:
        SAMTOOLS_VIEW_WRAPPER


rule cds_read_quant_tsv:
    input:
        PUFFERFISH_CDS_BAM_FILE,
    params:
        read_count=True,
    output:
        PUFFERFISH_CDS_READ_QUANT_TSV_FILE,
    log:
        PUFFERFISH_CDS_READ_QUANT_TSV_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_BAM_WRAPPER


rule cds_read_quant_hdf:
    conda:
        "../envs/vaex.yaml"
    input:
        PUFFERFISH_CDS_READ_QUANT_TSV_FILE,
    output:
        PUFFERFISH_CDS_READ_QUANT_HDF_FILE,
    log:
        PUFFERFISH_CDS_READ_QUANT_HDF_LOG,
    threads: CDS_READ_QUANT_THREADS
    script:
        "../scripts/tsv_to_hdf.py"


rule cedar_read_quant_tsv:
    input:
        pam=PUFFERFISH_GENOMIC_PAM_FILE,
        taxtree=NCBI_TAXDUMP_FILTERED_NODE_FILE,
        ref2tax=NCBI_ACC2TAXID_FILTERED_FILE,
    params:
        cedar=abspath(join(config["pufferfish"]["bin_dir"], "cedar")),
        extra=config["pufferfish"]["cedar"]["extra"],
    output:
        CEDAR_READ_QUANT_TSV_FILE,
    log:
        CEDAR_READ_QUANT_TSV_LOG,
    threads: CEDAR_READ_QUANT_THREADS
    wrapper:
        CEDAR_READ_QUANT_WRAPPER


rule cedar_read_quant_hdf:
    conda:
        "../envs/vaex.yaml"
    input:
        CEDAR_READ_QUANT_TSV_FILE,
    output:
        CEDAR_READ_QUANT_HDF_FILE,
    log:
        CEDAR_READ_QUANT_HDF_LOG,
    threads: CEDAR_READ_QUANT_THREADS
    script:
        "../scripts/tsv_to_hdf.py"
