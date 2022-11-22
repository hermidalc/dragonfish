from os.path import abspath


rule pufferfish_index:
    input:
        ref=REF_DEDUPED_ID_WITH_DECOY_FASTA_FILE,
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


rule diamond_uniprot_db:
    input:
        fastas=expand(UNIPROT_KB_FASTA_FILE, zip, **EXPAND_PARAMS),
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
