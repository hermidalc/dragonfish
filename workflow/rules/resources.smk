import pandas as pd


rule get_ncbi_acc2taxid_gz:
    params:
        NCBI_ACC2TAXID_GZ_URL,
    output:
        NCBI_ACC2TAXID_GZ_FILE,
    log:
        NCBI_ACC2TAXID_GZ_LOG,
    script:
        "../scripts/get_url_file.py"


rule gunzip_acc2taxid:
    input:
        NCBI_ACC2TAXID_GZ_FILE,
    output:
        NCBI_ACC2TAXID_FILE,
    log:
        NCBI_ACC2TAXID_LOG,
    script:
        "../scripts/unzip_file.py"


rule get_ncbi_taxdump_zip:
    params:
        NCBI_TAXDUMP_ZIP_URL,
    output:
        NCBI_TAXDUMP_ZIP_FILE,
    log:
        NCBI_TAXDUMP_ZIP_LOG,
    script:
        "../scripts/get_url_file.py"


rule unzip_ncbi_taxdump:
    input:
        NCBI_TAXDUMP_ZIP_FILE,
    output:
        NCBI_TAXDUMP_FILES,
    log:
        NCBI_TAXDUMP_FILES_LOG,
    script:
        "../scripts/unzip_file.py"


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
    input:
        NCBI_ASSEMBLY_SUMMARY_FILES,
    output:
        NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
    log:
        NCBI_MERGED_ASSEMBLY_SUMMARY_LOG,
    script:
        "../scripts/merge_ncbi_assembly_summaries.py"


rule get_uniprot_proteome_metadata:
    params:
        UNIPROT_PROTEOME_METADATA_URL,
    output:
        UNIPROT_PROTEOME_METADATA_FILE,
    log:
        UNIPROT_PROTEOME_METADATA_LOG,
    run:
        pd.read_csv(input[0], sep="\t").to_csv(output[0], sep="\t", index=False)


checkpoint get_ncbi_assembly_gz_files:
    conda:
        "../envs/get_ncbi_assembly_gz_files.yaml"
    input:
        proteome_file=UNIPROT_PROTEOME_METADATA_FILE,
        summary_file=NCBI_MERGED_ASSEMBLY_SUMMARY_FILE,
    params:
        file_exts=NCBI_ASSEMBLY_GZ_FILE_EXTS,
        parallel_backend=config["assembly"]["gz_file_job_backend"],
        verbosity=config["assembly"]["gz_file_job_verbosity"],
        debug=config["assembly"]["gz_file_job_debug"],
    output:
        directory(NCBI_ASSEMBLY_DIR),
    log:
        NCBI_ASSEMBLY_GZ_FILES_LOG,
    threads: 10
    script:
        "../scripts/get_ncbi_assembly_gz_files.py"


rule gunzip_ncbi_assembly_gz_file:
    input:
        NCBI_ASSEMBLY_GZ_FILE,
    output:
        NCBI_ASSEMBLY_FILE,
    log:
        NCBI_ASSEMBLY_FILE_LOG,
    script:
        "../scripts/unzip_file.py"


def gather_ncbi_assembly_files(wildcards):
    out_dir = checkpoints.get_ncbi_assembly_gz_files.get(**wildcards).output[0]
    basenames, exts = glob_wildcards(join(out_dir, "{asm_basename}.{asm_ext}.gz"))
    return expand(
        f"{out_dir}/{{asm_basename}}.{{asm_ext}}",
        zip,
        asm_basename=basenames,
        asm_ext=exts,
    )


rule finished:
    input:
        gather_ncbi_assembly_files,
    output:
        join(NCBI_GENOME_DIR, "aggregated.txt"),
    shell:
        """
        echo {input} > {output}
        """
