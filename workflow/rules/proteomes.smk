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


checkpoint split_uniprot_kb:
    conda:
        "../envs/lxml.yaml"
    input:
        UNIPROT_KB_FILE,
    params:
        basename="{ukb_basename}",
        split_size=lambda w: config["uniprot"]["kb"]["parse"]["split_size"][
            w.ukb_basename.split("_")[1]
        ],
        parser=config["uniprot"]["kb"]["parse"]["parser"],
    output:
        directory(UNIPROT_KB_SPLIT_DIR),
    log:
        UNIPROT_KB_SPLIT_LOG,
    script:
        "../scripts/split_uniprot_xml_file.py"


def gather_uniprot_kb_type_split_metadata_files(wildcards):
    split_dir = checkpoints.split_uniprot_kb.get(**wildcards).output[0]
    basenames, nums = glob_wildcards(
        join(split_dir, "{ukb_basename}_{ukb_snum}.xml.gz")
    )
    return sorted(
        expand(
            join(split_dir, "{ukb_basename}_{ukb_mtype}_{ukb_snum}.tsv"),
            zip,
            ukb_basename=basenames,
            ukb_snum=nums,
            allow_missing=True,
        )
    )


rule create_uniprot_kb_split_metadata:
    input:
        kb_file=UNIPROT_KB_SPLIT_FILE,
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
    params:
        dbxref_names=config["uniprot"]["kb"]["dbxref"]["names"],
    output:
        main=UNIPROT_KB_MAIN_METADATA_FILE,
        dbxref=UNIPROT_KB_DBXREF_METADATA_FILE,
    log:
        UNIPROT_KB_METADATA_LOG,
    script:
        "../scripts/create_uniprot_kb_metadata.py"


rule merge_uniprot_kb_type_metadata_splits:
    input:
        gather_uniprot_kb_type_split_metadata_files,
    output:
        UNIPROT_KB_TYPE_MERGED_METADATA_FILE,
    log:
        UNIPROT_KB_TYPE_MERGED_METADATA_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"


rule merge_uniprot_kb_metadata:
    input:
        sorted(
            expand(
                UNIPROT_KB_TYPE_MERGED_METADATA_FILE,
                ukb_basename=EXPAND_PARAMS["ukb_basename"],
                allow_missing=True,
            ),
            reverse=True,
        ),
    output:
        UNIPROT_KB_MERGED_METADATA_FILE,
    log:
        UNIPROT_KB_MERGED_METADATA_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"


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
