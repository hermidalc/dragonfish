from os.path import abspath


rule pufferfish_index:
    input:
        ref=REF_DEDUPED_ID_WITH_DECOY_FASTA_FILE,
        # decoys=GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE,
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


rule diamond_uniprot_db:
    input:
        fastas=UNIPROT_KB_MERGED_FASTA_FILE,
        taxonmap=UNIPROT_KB_TAXID_MAP_FILE,
        taxonnodes=NCBI_TAXDUMP_NODE_FILE,
        taxonnames=NCBI_TAXDUMP_NAME_FILE,
    output:
        DIAMOND_UNIPROT_DB_FILE,
    log:
        DIAMOND_UNIPROT_DB_LOG,
    threads: DIAMOND_THREADS
    wrapper:
        DIAMOND_MAKEDB_WRAPPER


rule paladin_uniprot_index:
    input:
        UNIPROT_KB_MERGED_FASTA_FILE,
    params:
        ref_type=3,
    output:
        directory(PALADIN_UNIPROT_INDEX_DIR),
    log:
        PALADIN_UNIPROT_INDEX_LOG,
    threads: PALADIN_THREADS
    wrapper:
        PALADIN_INDEX_WRAPPER
