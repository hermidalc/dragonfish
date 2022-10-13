rule get_ncbi_assembly_summary:
    params:
        NCBI_ASSEMBLY_SUMMARY_FILE_URL,
    output:
        NCBI_ASSEMBLY_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_SUMMARY_LOG,
    message:
        "{params}"
    retries: config["download"]["retries"]
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
        proteomes=UNIPROT_PROTEOME_METADATA_FILE,
        summary=NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
    params:
        genome_ext=config["ncbi"]["assembly"]["file"]["ext"]["genome"],
        cds_ext=config["ncbi"]["assembly"]["file"]["ext"]["cds"],
        skip=config["ncbi"]["assembly"]["file"]["download"]["skip"],
        md5_name=config["ncbi"]["assembly"]["file"]["download"]["md5_name"],
        retries=config["ncbi"]["assembly"]["file"]["download"]["file_retries"],
        retry_wait=config["ncbi"]["assembly"]["file"]["download"]["file_retry_wait"],
        backend=config["joblib"]["backend"],
        verbosity=config["joblib"]["verbosity"],
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_FILES_LOG,
    retries: config["download"]["retries"]
    threads: NCBI_ASSEMBLY_FILE_DOWNLOAD_THREADS
    script:
        "../scripts/get_ncbi_assemblies.py"


def gather_ncbi_assembly_fasta_files(wildcards):
    ncbi_assembly_dir = checkpoints.get_ncbi_assemblies.get(**wildcards).output[0]
    # XXX: workaround using separate unused asm_name wildcard for snakemake issue
    # #1849, using two of same wildcards in glob_wildcards pattern seems to be broken,
    # even though for regular rules it works. Also couldn't get asm_name wildcard
    # constraint regex patterns to work
    dirs, names = glob_wildcards(
        join(ncbi_assembly_dir, "{asm_dir}", "{asm_name}.fna.gz")
    )
    # required since some genome assemblies do not have cds_from_genomic files and
    # cannot use these dirs in expand for cds_from_genomic(_filtered) create fasta
    # list input, also cds_from_genomic_filtered files don't exist yet when this
    # function is called so need to use cds_from_genomic dirs for this call
    dirs = [
        d
        for d, n in zip(dirs, names)
        if n.replace(f"{d}_", "", 1)
        == (
            "cds_from_genomic"
            if wildcards.asm_type.startswith("cds_from_genomic")
            else wildcards.asm_type
        )
    ]
    return expand(
        f"{ncbi_assembly_dir}/{{asm_dir}}/{{asm_dir}}_{wildcards.asm_type}.fna.gz",
        asm_dir=dirs,
    )


rule create_ncbi_assembly_cds_filtered_fasta:
    input:
        NCBI_ASSEMBLY_CDS_FASTA_FILE,
    params:
        cmd="grep",
        pattern=config["ncbi"]["assembly"]["seqkit"]["grep"]["pattern"],
        extra=config["ncbi"]["assembly"]["seqkit"]["grep"]["extra"],
    output:
        NCBI_ASSEMBLY_CDS_FILTERED_FASTA_FILE,
    log:
        NCBI_ASSEMBLY_CDS_FILTERED_FASTA_LOG,
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
