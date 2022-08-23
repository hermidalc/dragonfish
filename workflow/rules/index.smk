rule create_pufferfish_reference_fasta_list:
    input:
        NCBI_ASSEMBLY_FASTA_LIST_FILE,
        GENCODE_GENOME_FIXED_FASTA_LIST_FILE,
    output:
        PUFFERFISH_REFERENCE_FASTA_LIST_FILE,
    log:
        PUFFERFISH_REFERENCE_FASTA_LIST_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"


rule create_pufferfish_reference_fasta:
    input:
        list_file=PUFFERFISH_REFERENCE_FASTA_LIST_FILE,
    params:
        extra=(
            " --only-id"
            f" --id-regexp '{PUFFERFISH_REFERENCE_SEQKIT_SEQ_ID_REGEX}'"
            f" --line-width {PUFFERFISH_REFERENCE_SEQKIT_SEQ_LINE_WIDTH}"
            " " + config["pufferfish"]["ref"]["seqkit"]["seq"]["extra_params"]
        ),
    output:
        PUFFERFISH_REFERENCE_FASTA_FILE,
    log:
        PUFFERFISH_REFERENCE_FASTA_LOG,
    threads: config["seqkit"]["seq"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER


rule create_pufferfish_index:
    input:
        ref=PUFFERFISH_REFERENCE_FASTA_FILE,
        decoys=GENCODE_GENOME_MERGED_FIXED_SEQ_ID_FILE,
    params:
        pufferfish=config["pufferfish"]["binary"],
        extra=config["pufferfish"]["index"]["extra_params"],
    output:
        directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    resources:
        tmpdir=TEMP_DIR,
    threads: PUFFERFISH_INDEX_THREADS
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
