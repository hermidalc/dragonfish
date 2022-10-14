from os.path import abspath


rule create_pufferfish_ref_fasta:
    input:
        list_file=NCBI_ASSEMBLY_FASTA_LIST_FILE,
    params:
        cmd="seq",
        extra=(
            "--only-id"
            + str(
                lambda wc: (
                    " --id-regexp "
                    + config["pufferfish"]["ref"]["seqkit"]["seq"][
                        "protein_id_pattern"
                    ]
                    if wc.asm_type.startswith("cds_from_genomic")
                    else ""
                )
            )
            + f' --line-width {config["seqkit"]["line_width"]}'
            f' {config["pufferfish"]["ref"]["seqkit"]["seq"]["extra"]}'
        ),
    output:
        PUFFERFISH_REF_FASTA_FILE,
    log:
        PUFFERFISH_REF_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule merge_pufferfish_ref_fastas:
    conda:
        "../envs/pigz.yaml"
    input:
        expand(PUFFERFISH_REF_FASTA_FILE, zip, **EXPAND_PARAMS)
        + expand(GENCODE_GENOME_FIXED_FASTA_FILE, zip, **EXPAND_PARAMS),
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


rule create_pufferfish_ref_merged_deduped_id_fasta:
    input:
        PUFFERFISH_REF_MERGED_FASTA_FILE,
    params:
        cmd="rename",
        extra=(
            f' --line-width {config["seqkit"]["line_width"]}'
            f' {config["pufferfish"]["ref"]["seqkit"]["rename"]["extra"]}'
        ),
    output:
        PUFFERFISH_REF_MERGED_DEDUPED_ID_FASTA_FILE,
    log:
        PUFFERFISH_REF_MERGED_DEDUPED_ID_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule create_pufferfish_index:
    input:
        ref=PUFFERFISH_REF_MERGED_DEDUPED_ID_FASTA_FILE,
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
