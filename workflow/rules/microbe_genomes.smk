rule ncbi_assembly_summary:
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
        "../scripts/url_file.py"


rule ncbi_assembly_merged_summary:
    conda:
        "../envs/pandas.yaml"
    input:
        NCBI_ASSEMBLY_SUMMARY_FILES,
    output:
        NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_MERGED_SUMMARY_LOG,
    script:
        "../scripts/ncbi_assembly_merged_summary.py"


rule ncbi_assembly_filtered_summary:
    conda:
        "../envs/pandas.yaml"
    input:
        summary=NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
        proteomes=UNIPROT_PROTEOMES_FILE,
    params:
        skip=config["ncbi"]["assembly"]["file"]["download"]["skip"],
    output:
        NCBI_ASSEMBLY_FILTERED_SUMMARY_FILE,
    log:
        NCBI_ASSEMBLY_FILTERED_SUMMARY_LOG,
    script:
        "../scripts/ncbi_assembly_filtered_summary.py"


checkpoint ncbi_assemblies:
    conda:
        "../envs/joblib.yaml"
    input:
        NCBI_ASSEMBLY_FILTERED_SUMMARY_FILE,
    params:
        file_exts=config["ncbi"]["assembly"]["file"]["exts"],
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
        "../scripts/ncbi_assemblies.py"


def gather_ncbi_assembly_fasta_files(wildcards):
    ncbi_assembly_dir = checkpoints.ncbi_assemblies.get(**wildcards).output[0]
    # XXX: workaround using separate unused asm_name wildcard for snakemake issue
    # #1849, using two of same wildcards in glob_wildcards pattern seems to be broken,
    # even though for regular rules it works. Also couldn't get asm_name wildcard
    # constraint regex patterns to work
    dirs, names = glob_wildcards(
        join(ncbi_assembly_dir, "{asm_dir}", "{asm_name}.fna.gz")
    )
    return sorted(
        expand(
            f"{ncbi_assembly_dir}/{{asm_dir}}/{{asm_dir}}_{wildcards.asm_type}.fna.gz",
            asm_dir=dirs,
        )
    )


rule ncbi_assembly_fasta_list:
    conda:
        "../envs/pandas.yaml"
    input:
        gather_ncbi_assembly_fasta_files,
    output:
        NCBI_ASSEMBLY_FASTA_LIST_FILE,
    log:
        NCBI_ASSEMBLY_FASTA_LIST_LOG,
    script:
        "../scripts/file_list_from_paths.py"


def gather_ncbi_assembly_cds_gtf_files(wildcards):
    ncbi_assembly_dir = checkpoints.ncbi_assemblies.get(**wildcards).output[0]
    # XXX: workaround here too
    dirs, names = glob_wildcards(
        join(ncbi_assembly_dir, "{asm_dir}", "{asm_name}_genomic.gff.gz")
    )
    return sorted(
        expand(
            f"{ncbi_assembly_dir}/{{asm_dir}}/{{asm_dir}}_genomic_cds.gtf",
            asm_dir=dirs,
        )
    )


rule ncbi_assembly_cds_gtf:
    conda:
        "../envs/gffutils.yaml"
    input:
        NCBI_ASSEMBLY_GFF_FILE,
    output:
        NCBI_ASSEMBLY_CDS_GTF_FILE,
    log:
        NCBI_ASSEMBLY_CDS_GTF_LOG,
    script:
        "../scripts/cds_gtf.py"


rule ncbi_assembly_merged_cds_gtf:
    input:
        gather_ncbi_assembly_cds_gtf_files,
    params:
        header="#gtf-version 2.2",
    output:
        NCBI_ASSEMBLY_MERGED_CDS_GTF_FILE,
    log:
        NCBI_ASSEMBLY_MERGED_CDS_GTF_LOG,
    script:
        "../scripts/merged_gff.py"


rule ncbi_assembly_taxid_map:
    conda:
        "../envs/pandas.yaml"
    input:
        files=expand(NCBI_ACC2TAXID_FILE, zip, **EXPAND_PARAMS),
        summary=NCBI_ASSEMBLY_FILTERED_SUMMARY_FILE,
    output:
        NCBI_ASSEMBLY_TAXID_MAP_FILE,
    log:
        NCBI_ASSEMBLY_TAXID_MAP_LOG,
    script:
        "../scripts/ncbi_assembly_taxid_map.py"
