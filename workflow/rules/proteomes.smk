rule uniprot_proteomes:
    params:
        tax_level=config["ncbi"]["taxonomy"]["level"],
        low_quality_pattern=config["ncbi"]["taxonomy"]["low_quality_pattern"],
        eukaryote_genera=config["ncbi"]["taxonomy"]["eukaryote_genera"],
        ref_url=config["uniprot"]["proteome"]["url"]["ref"],
        other_url=config["uniprot"]["proteome"]["url"]["other"],
        filter_domains=config["uniprot"]["proteome"]["filter_domains"],
        skip=config["uniprot"]["proteome"]["skip"],
        n_sample=config["uniprot"]["proteome"]["n_sample"],
        random_seed=config["random_seed"],
    output:
        UNIPROT_PROTEOMES_FILE,
    log:
        UNIPROT_PROTEOMES_LOG,
    retries: config["download"]["retries"]
    conda:
        "../envs/pandas.yaml"
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


rule uniprot_kb_merged_fasta:
    input:
        expand(UNIPROT_KB_FASTA_FILE, zip, **EXPAND_PARAMS),
    output:
        temp(UNIPROT_KB_MERGED_FASTA_FILE),
    log:
        UNIPROT_KB_MERGED_FASTA_LOG,
    # decompress takes ~1 thread in this context subtract 1
    threads: PIGZ_THREADS - 1
    conda:
        "../envs/pigz.yaml"
    shell:
        # creates a smaller gzip file than gzip cat
        # don't specify threads for decompress
        "pigz -dc {input} | pigz -p {threads} 1> {output} 2> {log}"


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
    input:
        kb=UNIPROT_KB_XML_FILE,
        proteomes=UNIPROT_PROTEOMES_FILE,
        split_pos=UNIPROT_KB_XML_SPLIT_POS_FILE,
    params:
        dbs=config["uniprot"]["kb"]["dbxref"]["dbs"],
        split_num="{ukb_snum}",
    output:
        temp(UNIPROT_KB_DBXREF_SPLIT_FILE),
    log:
        UNIPROT_KB_DBXREF_SPLIT_LOG,
    conda:
        "../envs/dbxref.yaml"
    script:
        "../scripts/uniprot_kb_dbxref_split.py"


rule uniprot_kb_merged_dbxref:
    input:
        UNIPROT_KB_DBXREF_SPLIT_FILES,
    output:
        UNIPROT_KB_DBXREF_FILE,
    log:
        UNIPROT_KB_DBXREF_LOG,
    conda:
        "../envs/vaex.yaml"
    script:
        "../scripts/merged_hdf.py"


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
    input:
        UNIPROT_KB_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["kb"]["idmap"]["parse"]["split_size"],
    output:
        UNIPROT_KB_GENBANK_IDMAP_FILE,
    log:
        UNIPROT_KB_GENBANK_IDMAP_LOG,
    conda:
        "../envs/vaex.yaml"
    script:
        "../scripts/uniprot_kb_genbank_idmap.py"


rule uniprot_kb_genbank_idmap_dbxrefs:
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
    conda:
        "../envs/vaex.yaml"
    script:
        "../scripts/uniprot_kb_genbank_idmap_dbxrefs.py"


rule uniprot_kb_taxid_map:
    input:
        UNIPROT_KB_IDMAP_FILE,
    output:
        UNIPROT_KB_TAXID_MAP_FILE,
    log:
        UNIPROT_KB_TAXID_MAP_LOG,
    script:
        "../scripts/uniprot_kb_taxid_map.py"
