rule get_uniprot_proteome_metadata:
    conda:
        "../envs/pandas.yaml"
    params:
        ref_url=config["uniprot"]["proteome"]["url"]["ref"],
        other_url=config["uniprot"]["proteome"]["url"]["other"],
        low_quality_pattern=config["ncbi"]["taxonomy"]["low_quality_pattern"],
        eukaryote_genera=config["ncbi"]["taxonomy"]["eukaryote_genera"],
        n_sample=config["uniprot"]["proteome"]["n_sample"],
        random_seed=config["uniprot"]["proteome"]["random_seed"],
        bacteria_only=config["uniprot"]["proteome"]["bacteria_only"],
    output:
        UNIPROT_PROTEOME_METADATA_FILE,
    log:
        UNIPROT_PROTEOME_METADATA_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/get_uniprot_proteome_metadata.py"


rule get_uniprot_kb:
    params:
        UNIPROT_KB_FILE_URL,
    output:
        UNIPROT_KB_FILE,
    log:
        UNIPROT_KB_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/get_url_file.py"


rule get_uniprot_kb_split_pos:
    input:
        UNIPROT_KB_FILE,
    params:
        kb_size=lambda wc: config["uniprot"]["kb"]["kb_sizes"][
            EXPAND_PARAMS["ukb_basename"].index(wc.ukb_basename)
        ],
        split_size=lambda wc: config["uniprot"]["kb"]["parse"]["split_sizes"][
            EXPAND_PARAMS["ukb_basename"].index(wc.ukb_basename)
        ],
    output:
        UNIPROT_KB_SPLIT_POS_FILE,
    log:
        UNIPROT_KB_SPLIT_POS_LOG,
    script:
        "../scripts/get_uniprot_kb_split_pos.py"


rule create_uniprot_kb_dbxref:
    conda:
        "../envs/dbxref.yaml"
    input:
        kb=UNIPROT_KB_FILE,
        proteomes=UNIPROT_PROTEOME_METADATA_FILE,
        split_pos=UNIPROT_KB_SPLIT_POS_FILE,
    params:
        dbs=config["uniprot"]["kb"]["dbxref"]["dbs"],
        split_num="{ukb_snum}",
    output:
        temp(UNIPROT_KB_DBXREF_FILE),
    log:
        UNIPROT_KB_DBXREF_LOG,
    script:
        "../scripts/create_uniprot_kb_dbxref.py"


rule merge_uniprot_kb_dbxrefs:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_DBXREF_FILES,
    output:
        UNIPROT_KB_DBXREF_HDF_FILE,
    log:
        UNIPROT_KB_DBXREF_HDF_LOG,
    script:
        "../scripts/merge_uniprot_kb_dbxrefs.py"


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


rule create_uniprot_kb_idmap:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["kb"]["idmap"]["parse"]["split_size"],
    output:
        UNIPROT_KB_IDMAP_HDF_FILE,
    log:
        UNIPROT_KB_IDMAP_HDF_LOG,
    resources:
        tmpdir=TEMP_DIR,
    script:
        "../scripts/create_uniprot_kb_idmap.py"


rule create_uniprot_kb_idmap_dbxref:
    conda:
        "../envs/vaex.yaml"
    input:
        idmap=UNIPROT_KB_IDMAP_HDF_FILE,
        dbxref=UNIPROT_KB_DBXREF_HDF_FILE,
    params:
        db=lambda wc: config["uniprot"]["kb"]["dbxref"]["dbs"][
            EXPAND_PARAMS["ukb_dbxref_db"].index(wc.ukb_dbxref_db)
        ],
    output:
        UNIPROT_KB_IDMAP_DBXREF_HDF_FILE,
    log:
        UNIPROT_KB_IDMAP_DBXREF_HDF_LOG,
    threads: UNIPROT_KB_IDMAP_DBXREF_THREADS
    script:
        "../scripts/create_uniprot_kb_idmap_dbxref.py"
