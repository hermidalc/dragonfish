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
    output:
        UNIPROT_PROTEOME_METADATA_FILE,
    log:
        UNIPROT_PROTEOME_METADATA_LOG,
    script:
        "../scripts/get_uniprot_proteome_metadata.py"


rule get_uniprot_kb_file:
    params:
        UNIPROT_KB_FILE_URL,
    output:
        UNIPROT_KB_FILE,
    log:
        UNIPROT_KB_LOG,
    message:
        "{params}"
    script:
        "../scripts/get_url_file.py"
