rule ref_fasta:
    input:
        list_file=NCBI_ASSEMBLY_FASTA_LIST_FILE,
    params:
        extra="--only-id " + config["ref"]["seqkit"]["seq"]["extra"],
    output:
        REF_FASTA_FILE,
    log:
        REF_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER


rule ref_deduped_id_fasta:
    input:
        REF_FASTA_FILE,
    params:
        extra=config["ref"]["seqkit"]["rename"]["extra"],
    output:
        REF_DEDUPED_ID_FASTA_FILE,
    log:
        REF_DEDUPED_ID_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_RENAME_WRAPPER


rule ref_deduped_id_with_decoys_fasta:
    input:
        ref=REF_DEDUPED_ID_FASTA_FILE,
        decoys=[GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE, T2T_GENOME_FASTA_FILE],
    output:
        REF_DEDUPED_ID_WITH_DECOY_FASTA_FILE,
    log:
        REF_DEDUPED_ID_WITH_DECOY_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: PIGZ_THREADS - 1
    conda:
        "../envs/pigz.yaml"
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input.ref} {input.decoys} | pigz -p {threads} 1> {output} 2> {log}"
