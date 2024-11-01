rule host_genome_index:
    input:
        fasta=T2T_GENOME_FASTA_FILE,
    params:
        prefix=HOST_GENOME_INDEX_PREFIX,
        extra=f"{config[HOST_FILTER_MODE]['build']['extra']} --seed {config['random_seed']}",
    output:
        directory(HOST_GENOME_INDEX_DIR),
    log:
        HOST_GENOME_INDEX_LOG,
    threads: HOST_BUILD_THREADS
    wrapper:
        HOST_BUILD_WRAPPER


# rule host_filtered_fastq_pe:
#     input:
#         reads=[GDC_UNMAPPED_FASTQ_R1_FILE, GDC_UNMAPPED_FASTQ_R2_FILE],
#         dir=HOST_GENOME_INDEX_DIR,
#     params:
#         idx=HOST_GENOME_INDEX_PREFIX,
#         extra=f"{config[HOST_FILTER_MODE]['align']['extra']} --seed {config['random_seed']}",
#     output:
#         # output=temp(HOST_BAM_PE_FILE),
#         unconcordant=[
#             temp(HOST_FILTERED_FASTQ_R1_FILE),
#             temp(HOST_FILTERED_FASTQ_R2_FILE),
#         ],
#     log:
#         HOST_FILTERED_FASTQ_LOG,
#     threads: HOST_ALIGN_THREADS
#     wrapper:
#         HOST_ALIGN_WRAPPER
# rule host_filtered_fastq_se:
#     input:
#         reads=[GDC_UNMAPPED_FASTQ_SE_FILE],
#         dir=HOST_GENOME_INDEX_DIR,
#     params:
#         idx=HOST_GENOME_INDEX_PREFIX,
#         extra=f"{config[HOST_FILTER_MODE]['align']['extra']} --seed {config['random_seed']}",
#     output:
#         # output=temp(HOST_BAM_SE_FILE),
#         unaligned=temp(HOST_FILTERED_FASTQ_SE_FILE),
#     log:
#         HOST_FILTERED_FASTQ_LOG,
#     threads: HOST_ALIGN_THREADS
#     wrapper:
#         HOST_ALIGN_WRAPPER
