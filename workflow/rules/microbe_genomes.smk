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


checkpoint get_ncbi_assembly_files:
    conda:
        "../envs/joblib.yaml"
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=NCBI_ASSEMBLY_MERGED_SUMMARY_FILE,
    params:
        file_exts=NCBI_ASSEMBLY_FILE_EXTS,
        backend=config["ncbi"]["assembly"]["download_backend"],
        verbosity=config["ncbi"]["assembly"]["download_verbosity"],
        debug=config["ncbi"]["assembly"]["download_debug"],
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_FILES_LOG,
    threads: config["ncbi"]["assembly"]["download_threads"]
    script:
        "../scripts/get_ncbi_assembly_files.py"


def gather_ncbi_assembly_fasta_files(wildcards):
    out_dir = checkpoints.get_ncbi_assembly_files.get(**wildcards).output[0]
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


rule create_ncbi_assembly_fasta_list_file:
    conda:
        "../envs/pandas.yaml"
    input:
        gather_ncbi_assembly_fasta_files,
    output:
        NCBI_ASSEMBLY_FASTA_LIST_FILE,
    script:
        "../scripts/create_ncbi_assembly_fasta_list_file.py"


rule create_ncbi_reference_fasta:
    input:
        list_file=NCBI_ASSEMBLY_FASTA_LIST_FILE,
    params:
        extra=(
            " --only-id"
            f" --id-regexp '{NCBI_ASSEMBLY_FASTA_SEQKIT_SEQ_ID_REGEX}'"
            f" --line-width {SEQKIT_FASTA_LINE_WIDTH}"
            " " + config["seqkit"]["seq"]["extra_params"]
        ),
    output:
        NCBI_REFERENCE_FASTA_FILE,
    log:
        NCBI_REFERENCE_FASTA_LOG,
    threads: config["seqkit"]["seq"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER
