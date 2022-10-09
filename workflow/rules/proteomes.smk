rule get_uniprot_proteome_metadata:
    conda:
        "../envs/pandas.yaml"
    params:
        ref_url=config["uniprot"]["proteome"]["url"]["ref"],
        other_url=config["uniprot"]["proteome"]["url"]["other"],
        low_quality_regex=config["ncbi"]["taxonomy"]["low_quality_regex"],
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


rule get_uniprot_kb_dbxrefs:
    conda:
        "../envs/biopython.yaml"
    input:
        kb=UNIPROT_KB_FILE,
        proteome=UNIPROT_PROTEOME_METADATA_FILE,
        split_pos=UNIPROT_KB_SPLIT_POS_FILE,
    params:
        dbxref_names=config["uniprot"]["kb"]["dbxref"]["names"],
        split_num="{ukb_snum}",
    output:
        UNIPROT_KB_DBXREF_FILE,
    log:
        UNIPROT_KB_DBXREF_LOG,
    script:
        "../scripts/get_uniprot_kb_dbxrefs.py"


rule create_uniprot_kb_dbxref_hdf:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_DBXREF_FILE,
    output:
        UNIPROT_KB_DBXREF_HDF_FILE,
    log:
        UNIPROT_KB_DBXREF_HDF_LOG,
    script:
        "../scripts/create_uniprot_kb_dbxref_hdf.py"


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


checkpoint create_uniprot_kb_idmap_hdf:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_KB_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["kb"]["idmap"]["parse"]["split_size"],
    output:
        directory(UNIPROT_KB_IDMAP_HDF_DIR),
    log:
        UNIPROT_KB_IDMAP_HDF_LOG,
    script:
        "../scripts/create_uniprot_kb_idmap_hdf.py"


def gather_uniprot_kb_idmap_hdf_files(wildcards):
    hdf_dir = checkpoints.create_uniprot_kb_idmap_hdf.get(**wildcards).output[0]
    file_wc_path = join(hdf_dir, f"{UNIPROT_KB_IDMAP_FILE_BASENAME}_{{i}}.hdf5")
    return sorted(expand(file_wc_path, i=glob_wildcards(file_wc_path).i))
