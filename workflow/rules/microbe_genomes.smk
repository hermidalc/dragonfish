rule get_ncbi_assembly_summary:
    params:
        NCBI_ASSEMBLY_SUMMARY_FILE_URL,
    output:
        NCBI_ASSEMBLY_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_SUMMARY_LOG,
    script:
        "../scripts/get_url_file.py"


rule merge_ncbi_assembly_summaries:
    conda:
        "../envs/pandas.yaml"
    input:
        NCBI_ASSEMBLY_SUMMARY_FILES,
    output:
        NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_MERGED_SUMMARY_LOG,
    script:
        "../scripts/merge_ncbi_assembly_summaries.py"


checkpoint get_ncbi_assemblies:
    conda:
        "../envs/joblib.yaml"
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
    params:
        genome_ext=config["ncbi"]["assembly"]["file"]["ext"]["genome"],
        cds_ext=config["ncbi"]["assembly"]["file"]["ext"]["cds"],
        skip=config["ncbi"]["assembly"]["file"]["download"]["skip"],
        md5_name=config["ncbi"]["assembly"]["file"]["download"]["md5_name"],
        retries=config["ncbi"]["assembly"]["file"]["download"]["file_retries"],
        retry_wait=config["ncbi"]["assembly"]["file"]["download"]["file_retry_wait"],
        backend=config["ncbi"]["assembly"]["file"]["download"]["backend"],
        verbosity=config["ncbi"]["assembly"]["file"]["download"]["verbosity"],
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_FILES_LOG,
    retries: config["ncbi"]["assembly"]["file"]["download"]["job_retries"]
    threads: NCBI_ASSEMBLY_FILE_DOWNLOAD_THREADS
    script:
        "../scripts/get_ncbi_assemblies.py"


def gather_ncbi_assembly_fasta_files(wildcards):
    ncbi_assembly_dir = checkpoints.get_ncbi_assemblies.get(**wildcards).output[0]
    # XXX: workaround for snakemake issue #1849, using two of same wildcard names in
    # glob_wildcards pattern doesn't work even though for regular rules it does and
    # also couldn't get negative assertion wildcard contraints on asm_dir to work
    dirs, _ = glob_wildcards(join(ncbi_assembly_dir, "{asm_dir}", "{asm_name}.fna.gz"))
    return expand(
        f"{ncbi_assembly_dir}/{{asm_dir}}/{{asm_dir}}_{wildcards.asm_type}.fna.gz",
        asm_dir=dirs,
    )


rule create_ncbi_assembly_cds_no_pseudo_fasta:
    input:
        NCBI_ASSEMBLY_CDS_FASTA_FILE,
    params:
        cmd="grep",
        pattern="\[pseudo=true\]",
        extra="--by-name --use-regexp --invert-match --ignore-case",
    output:
        NCBI_ASSEMBLY_CDS_NO_PSEUDO_FASTA_FILE,
    log:
        NCBI_ASSEMBLY_CDS_NO_PSEUDO_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule create_ncbi_assembly_fasta_list:
    conda:
        "../envs/pandas.yaml"
    input:
        gather_ncbi_assembly_fasta_files,
    params:
        sort_by=["dir", "name"],
    output:
        NCBI_ASSEMBLY_FASTA_LIST_FILE,
    log:
        NCBI_ASSEMBLY_FASTA_LIST_LOG,
    script:
        "../scripts/create_file_list_from_paths.py"
