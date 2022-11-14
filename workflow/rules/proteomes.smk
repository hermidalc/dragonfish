rule uniprot_proteomes:
    conda:
        "../envs/pandas.yaml"
    params:
        tax_level=config["ncbi"]["taxonomy"]["level"],
        low_quality_pattern=config["ncbi"]["taxonomy"]["low_quality_pattern"],
        eukaryote_genera=config["ncbi"]["taxonomy"]["eukaryote_genera"],
        ref_url=config["uniprot"]["proteome"]["url"]["ref"],
        other_url=config["uniprot"]["proteome"]["url"]["other"],
        filter_domains=config["uniprot"]["proteome"]["filter_domains"],
        n_sample=config["uniprot"]["proteome"]["n_sample"],
        random_seed=config["random_seed"],
    output:
        UNIPROT_PROTEOMES_FILE,
    log:
        UNIPROT_PROTEOMES_LOG,
    retries: config["download"]["retries"]
    script:
        "../scripts/uniprot_proteomes.py"


rule uniprot_kb_xml:
    params:
        UNIPROT_KB_XML_URL,
    output:
        UNIPROT_KB_XML_FILE,
    log:
        UNIPROT_KB_XML_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/url_file.py"


rule uniprot_kb_fasta:
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
        "../scripts/url_file.py"


rule uniprot_kb_xml_split_pos:
    input:
        UNIPROT_KB_XML_FILE,
    params:
        num_splits=lambda wc: config["uniprot"]["kb"]["parse"]["num_splits"][
            EXPAND_PARAMS["ukb_basename"].index(wc.ukb_basename)
        ],
        split_size=lambda wc: config["uniprot"]["kb"]["parse"]["split_size"],
    output:
        UNIPROT_KB_XML_SPLIT_POS_FILE,
    log:
        UNIPROT_KB_XML_SPLIT_POS_LOG,
    script:
        "../scripts/uniprot_kb_xml_split_pos.py"


rule uniprot_kb_dbxref_split:
    conda:
        "../envs/dbxref.yaml"
    input:
        kb=UNIPROT_KB_XML_FILE,
        proteomes=UNIPROT_PROTEOMES_FILE,
        split_pos=UNIPROT_KB_XML_SPLIT_POS_FILE,
    params:
        dbs=config["uniprot"]["kb"]["dbxref"]["dbs"],
        split_num="{ukb_snum}",
    output:
        UNIPROT_KB_DBXREF_SPLIT_FILE,
    log:
        UNIPROT_KB_DBXREF_SPLIT_LOG,
    script:
        "../scripts/uniprot_kb_dbxref_split.py"


rule uniprot_kb_merged_dbxref:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_DBXREF_SPLIT_FILES,
    output:
        UNIPROT_KB_DBXREF_FILE,
    log:
        UNIPROT_KB_DBXREF_LOG,
    script:
        "../scripts/merged_hdf_splits.py"


rule uniprot_kb_idmap:
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
        "../scripts/url_file.py"


rule uniprot_kb_genbank_idmap:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["kb"]["idmap"]["parse"]["split_size"],
    output:
        UNIPROT_KB_GENBANK_IDMAP_FILE,
    log:
        UNIPROT_KB_GENBANK_IDMAP_LOG,
    script:
        "../scripts/uniprot_kb_genbank_idmap.py"


rule uniprot_kb_genbank_idmap_dbxrefs:
    conda:
        "../envs/vaex.yaml"
    input:
        idmap=UNIPROT_KB_GENBANK_IDMAP_FILE,
        dbxref=UNIPROT_KB_DBXREF_FILE,
    params:
        db=lambda wc: config["uniprot"]["kb"]["dbxref"]["dbs"][
            EXPAND_PARAMS["ukb_dbxref_db"].index(wc.ukb_dbxref_db)
        ],
    output:
        UNIPROT_KB_GENBANK_IDMAP_DBXREF_FILE,
    log:
        UNIPROT_KB_GENBANK_IDMAP_DBXREF_LOG,
    threads: UNIPROT_KB_GENBANK_IDMAP_DBXREF_THREADS
    script:
        "../scripts/uniprot_kb_genbank_idmap_dbxrefs.py"
