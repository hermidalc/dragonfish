from os.path import abspath


rule pufferfish_ref_fasta:
    input:
        list_file=NCBI_ASSEMBLY_FASTA_LIST_FILE,
    params:
        id_regexp=lambda wc: (
            config["pufferfish"]["ref"]["seqkit"]["seq"]["protein_id_pattern"]
            if wc.asm_type.startswith("cds_from_genomic")
            else None
        ),
        extra="--only-id " + config["pufferfish"]["ref"]["seqkit"]["seq"]["extra"],
    output:
        PUFFERFISH_REF_BASE_FASTA_FILE,
    log:
        PUFFERFISH_REF_BASE_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER


rule pufferfish_ref_merged_fasta:
    conda:
        "../envs/pigz.yaml"
    input:
        expand(PUFFERFISH_REF_BASE_FASTA_FILE, zip, **EXPAND_PARAMS),
    output:
        PUFFERFISH_REF_MERGED_FASTA_FILE,
    log:
        PUFFERFISH_REF_MERGED_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: PUFFERFISH_REF_PIGZ_THREADS - 1
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input} | pigz -p {threads} 1> {output} 2> {log}"


rule pufferfish_ref_merged_decoy_fasta:
    conda:
        "../envs/pigz.yaml"
    input:
        [PUFFERFISH_REF_FASTA_FILE]
        + expand(GENCODE_GENOME_FIXED_FASTA_FILE, zip, **EXPAND_PARAMS),
    output:
        PUFFERFISH_REF_MERGED_DECOY_FASTA_FILE,
    log:
        PUFFERFISH_REF_MERGED_DECOY_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: PUFFERFISH_REF_PIGZ_THREADS - 1
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input} | pigz -p {threads} 1> {output} 2> {log}"


rule pufferfish_ref_deduped_id_fasta:
    input:
        PUFFERFISH_REF_MERGED_DECOY_FASTA_FILE,
    params:
        extra=config["pufferfish"]["ref"]["seqkit"]["rename"]["extra"],
    output:
        PUFFERFISH_REF_MERGED_DECOY_DEDUPED_ID_FASTA_FILE,
    log:
        PUFFERFISH_REF_MERGED_DECOY_DEDUPED_ID_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_RENAME_WRAPPER


rule pufferfish_index:
    input:
        ref=PUFFERFISH_REF_MERGED_DECOY_DEDUPED_ID_FASTA_FILE,
        decoys=GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE,
    params:
        pufferfish=abspath(join(config["pufferfish"]["bin_dir"], "pufferfish")),
        extra=config["pufferfish"]["index"]["extra"],
    output:
        directory(PUFFERFISH_INDEX_DIR),
    log:
        PUFFERFISH_INDEX_LOG,
    threads: PUFFERFISH_INDEX_THREADS
    wrapper:
        PUFFERFISH_INDEX_WRAPPER
