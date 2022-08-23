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
        genome_ext=config["ncbi"]["assembly"]["file"]["genome_ext"],
        other_exts=config["ncbi"]["assembly"]["file"]["other_exts"],
        skip=config["ncbi"]["assembly"]["file"]["download"]["skip"],
        md5_name=config["ncbi"]["assembly"]["file"]["download"]["md5_name"],
        retries=config["ncbi"]["assembly"]["file"]["download"]["retries"],
        force=config["ncbi"]["assembly"]["file"]["download"]["force"],
        backend=config["ncbi"]["assembly"]["file"]["download"]["backend"],
        verbosity=config["ncbi"]["assembly"]["file"]["download"]["verbosity"],
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_FILES_LOG,
    threads: NCBI_ASSEMBLY_FILE_DOWNLOAD_THREADS
    script:
        "../scripts/get_ncbi_assemblies.py"


def gather_ncbi_assembly_fasta_files(wildcards):
    out_dir = checkpoints.get_ncbi_assemblies.get(**wildcards).output[0]
    dirs, basenames, exts = glob_wildcards(
        join(out_dir, "{asm_dir}", "{asm_basename}.{asm_ext}.gz")
    )
    return expand(
        f"{out_dir}/{{asm_dir}}/{{asm_basename}}.{{asm_ext}}.gz",
        zip,
        asm_dir=dirs,
        asm_basename=basenames,
        asm_ext=exts,
    )


rule create_ncbi_assembly_fasta_list:
    conda:
        "../envs/pandas.yaml"
    input:
        gather_ncbi_assembly_fasta_files,
    params:
        sort_by=["dir", "name"],
        sort_ascending=[True, False],
    output:
        NCBI_ASSEMBLY_FASTA_LIST_FILE,
    log:
        NCBI_ASSEMBLY_FASTA_LIST_LOG,
    script:
        "../scripts/create_file_list_from_paths.py"
