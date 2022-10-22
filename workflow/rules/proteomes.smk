rule get_uniprot_proteomes:
    conda:
        "../envs/pandas.yaml"
    params:
        ref_url=config["uniprot"]["proteome"]["url"]["ref"],
        other_url=config["uniprot"]["proteome"]["url"]["other"],
        low_quality_pattern=config["ncbi"]["taxonomy"]["low_quality_pattern"],
        eukaryote_genera=config["ncbi"]["taxonomy"]["eukaryote_genera"],
        n_sample=config["uniprot"]["proteome"]["n_sample"],
        random_seed=config["uniprot"]["proteome"]["random_seed"],
        filter_domains=config["uniprot"]["proteome"]["filter_domains"],
    output:
        UNIPROT_PROTEOMES_FILE,
    log:
        UNIPROT_PROTEOMES_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/get_uniprot_proteomes.py"


rule get_uniprot_kb:
    params:
        UNIPROT_KB_URL,
    output:
        UNIPROT_KB_FILE,
    log:
        UNIPROT_KB_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/get_url_file.py"


rule get_uniprot_kb_fasta:
    params:
        UNIPROT_KB_FASTA_URL,
    output:
        UNIPROT_KB_FASTA_FILE,
    log:
        UNIPROT_KB_FASTA_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/get_url_file.py"


rule get_uniprot_kb_split_pos:
    input:
        UNIPROT_KB_FILE,
    params:
        num_splits=lambda wc: config["uniprot"]["kb"]["parse"]["num_splits"][
            EXPAND_PARAMS["ukb_basename"].index(wc.ukb_basename)
        ],
        split_size=lambda wc: config["uniprot"]["kb"]["parse"]["split_size"],
    output:
        UNIPROT_KB_SPLIT_POS_FILE,
    log:
        UNIPROT_KB_SPLIT_POS_LOG,
    script:
        "../scripts/get_uniprot_kb_split_pos.py"


rule create_uniprot_kb_dbxref_split:
    conda:
        "../envs/dbxref.yaml"
    input:
        kb=UNIPROT_KB_FILE,
        proteomes=UNIPROT_PROTEOMES_FILE,
        split_pos=UNIPROT_KB_SPLIT_POS_FILE,
    params:
        dbs=config["uniprot"]["kb"]["dbxref"]["dbs"],
        split_num="{ukb_snum}",
    output:
        UNIPROT_KB_DBXREF_SPLIT_FILE,
    log:
        UNIPROT_KB_DBXREF_SPLIT_LOG,
    script:
        "../scripts/create_uniprot_kb_dbxref_split.py"


rule merge_uniprot_kb_dbxrefs:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_DBXREF_SPLIT_FILES,
    output:
        UNIPROT_KB_DBXREF_FILE,
    log:
        UNIPROT_KB_DBXREF_LOG,
    script:
        "../scripts/merge_uniprot_kb_dbxref_splits.py"


rule get_uniprot_kb_idmap:
    params:
        UNIPROT_KB_IDMAP_FILE_URL,
    output:
        UNIPROT_KB_IDMAP_FILE,
    log:
        UNIPROT_KB_IDMAP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/get_url_file.py"


rule create_uniprot_kb_genbank_cds_idmap:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["kb"]["idmap"]["parse"]["split_size"],
    output:
        UNIPROT_KB_GENBANK_CDS_IDMAP_FILE,
    log:
        UNIPROT_KB_GENBANK_CDS_IDMAP_LOG,
    resources:
        tmpdir=TEMP_DIR,
    script:
        "../scripts/create_uniprot_kb_genbank_cds_idmap.py"


rule create_uniprot_kb_genbank_cds_idmap_dbxref:
    conda:
        "../envs/vaex.yaml"
    input:
        idmap=UNIPROT_KB_GENBANK_CDS_IDMAP_FILE,
        dbxref=UNIPROT_KB_DBXREF_FILE,
    params:
        db=lambda wc: config["uniprot"]["kb"]["dbxref"]["dbs"][
            EXPAND_PARAMS["ukb_dbxref_db"].index(wc.ukb_dbxref_db)
        ],
    output:
        UNIPROT_KB_GENBANK_CDS_IDMAP_DBXREF_FILE,
    log:
        UNIPROT_KB_GENBANK_CDS_IDMAP_DBXREF_LOG,
    threads: UNIPROT_KB_GENBANK_CDS_IDMAP_DBXREF_THREADS
    script:
        "../scripts/create_uniprot_kb_genbank_cds_idmap_dbxref.py"
