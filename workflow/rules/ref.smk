rule ref_assembly_fasta:
    input:
        list_file=NCBI_ASSEMBLY_FASTA_LIST_FILE,
    params:
        id_regexp=lambda wc: (
            config["ref"]["seqkit"]["seq"]["pattern"]
            if wc.asm_type.startswith("cds_from_genomic")
            else None
        ),
        extra="--only-id " + config["ref"]["seqkit"]["seq"]["extra"],
    output:
        REF_ASSEMBLY_FASTA_FILE,
    log:
        REF_ASSEMBLY_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER


rule ref_assembly_deduped_id_fasta:
    input:
        REF_ASSEMBLY_FASTA_FILE,
    params:
        extra=config["ref"]["seqkit"]["rename"]["extra"],
    output:
        REF_ASSEMBLY_DEDUPED_ID_FASTA_FILE,
    log:
        REF_ASSEMBLY_DEDUPED_ID_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_RENAME_WRAPPER


rule ref_assembly_deduped_qnames:
    input:
        REF_ASSEMBLY_DEDUPED_ID_FASTA_FILE,
    output:
        REF_ASSEMBLY_DEDUPED_QNAME_FILE,
    log:
        REF_ASSEMBLY_DEDUPED_QNAME_LOG,
    shell:
        "pigz -dc {input} | grep '^>' | cut -c2- | awk '{{print \"SN:\"$0}}' 1> {output} 2> {log}"


rule ref_merged_deduped_id_fasta:
    conda:
        "../envs/pigz.yaml"
    input:
        expand(REF_ASSEMBLY_DEDUPED_ID_FASTA_FILE, zip, **EXPAND_PARAMS),
    output:
        REF_MERGED_DEDUPED_ID_FASTA_FILE,
    log:
        REF_MERGED_DEDUPED_ID_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: REF_PIGZ_THREADS - 1
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input} | pigz -p {threads} 1> {output} 2> {log}"


rule ref_deduped_id_with_decoy_fasta:
    conda:
        "../envs/pigz.yaml"
    input:
        [REF_DEDUPED_ID_FASTA_FILE]
        + expand(GENCODE_GENOME_FIXED_FASTA_FILE, zip, **EXPAND_PARAMS),
    output:
        REF_DEDUPED_ID_WITH_DECOY_FASTA_FILE,
    log:
        REF_DEDUPED_ID_WITH_DECOY_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: REF_PIGZ_THREADS - 1
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input} | pigz -p {threads} 1> {output} 2> {log}"
