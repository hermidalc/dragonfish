rule get_uniprot_proteome_metadata:
    conda:
        "../envs/pandas.yaml"
    params:
        ref_url=config["uniprot"]["proteome"]["ref_metadata_url"],
        other_url=config["uniprot"]["proteome"]["other_metadata_url"],
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


rule create_uniprot_kb_dbxref_split:
    conda:
        "../envs/biopython.yaml"
    input:
        kb_file=UNIPROT_KB_FILE,
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
    params:
        dbxref_names=config["uniprot"]["kb"]["dbxref"]["names"],
        kb_size=lambda wc: config["uniprot"]["kb"]["kb_sizes"][
            0 if wc.ukb_basename == "uniprot_sprot" else 1
        ],
        split_size=config["uniprot"]["kb"]["parse"]["split_size"],
        split_num="{ukb_snum}",
    output:
        UNIPROT_KB_DBXREF_SPLIT_FILE,
    log:
        UNIPROT_KB_DBXREF_SPLIT_LOG,
    script:
        "../scripts/create_uniprot_kb_metadata.py"


rule merge_uniprot_kb_dbxref_splits:
    conda:
        "../envs/pigz.yaml"
    input:
        UNIPROT_KB_DBXREF_SPLIT_FILES,
    output:
        UNIPROT_KB_MERGED_DBXREF_FILE,
    log:
        UNIPROT_KB_MERGED_DBXREF_LOG,
    shell:
        "pigz -dc {input} | pigz -p 1 1> {output} 2> {log}"


rule get_uniprot_idmap:
    params:
        UNIPROT_IDMAP_FILE_URL,
    output:
        UNIPROT_IDMAP_FILE,
    log:
        UNIPROT_IDMAP_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
    script:
        "../scripts/get_url_file.py"


checkpoint create_uniprot_idmap_hdf:
    conda:
        "../envs/vaex.yaml"
    input:
        UNIPROT_IDMAP_FILE,
    params:
        split_size=config["uniprot"]["idmap"]["parse"]["split_size"],
    output:
        directory(UNIPROT_IDMAP_HDF_DIR),
    log:
        UNIPROT_IDMAP_HDF_LOG,
    script:
        "../scripts/create_uniprot_idmap_hdf.py"


def gather_uniprot_idmap_hdf_files(wildcards):
    hdf_dir = checkpoints.create_uniprot_idmap_hdf.get(**wildcards).output[0]
    file_wc_path = join(hdf_dir, f"{UNIPROT_IDMAP_FILE_BASENAME}_{{i}}.hdf5")
    return expand(file_wc_path, i=glob_wildcards(file_wc_path).i)
