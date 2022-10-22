rule gencode_genome_fasta:
    params:
        protocol=GENCODE_PROTOCOL,
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
        regions=GENCODE_REGIONS,
    output:
        GENCODE_GENOME_FASTA_FILE,
    resources:
        gencode_download_jobs=1,
    retries: config["download"]["retries"]
    wrapper:
        GENCODE_GENOME_SEQ_WRAPPER


rule gencode_genome_fasta_id_fix:
    input:
        GENCODE_GENOME_FASTA_FILE,
    params:
        cmd="replace",
        pattern=config["gencode"]["seqkit"]["replace"]["pattern"],
        replacement=lambda wc: (
            f"{wc.gc_species.title()}_{wc.gc_release}_{wc.gc_build}_${{1}}"
        ),
        extra=config["gencode"]["seqkit"]["replace"]["extra"],
    output:
        GENCODE_GENOME_FIXED_FASTA_FILE,
    log:
        GENCODE_GENOME_FIXED_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule gencode_genome_fixed_fasta_ids:
    input:
        GENCODE_GENOME_FIXED_FASTA_FILE,
    params:
        cmd="seq",
        extra="--name",
    output:
        GENCODE_GENOME_FIXED_FASTA_ID_FILE,
    log:
        GENCODE_GENOME_FIXED_FASTA_ID_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule gencode_genome_merged_fixed_fasta_ids:
    input:
        expand(GENCODE_GENOME_FIXED_FASTA_ID_FILE, zip, **EXPAND_PARAMS),
    output:
        GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE,
    log:
        GENCODE_GENOME_MERGED_FIXED_FASTA_ID_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"


rule gencode_genome_annot:
    params:
        protocol=GENCODE_PROTOCOL,
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        regions=GENCODE_REGIONS,
        annot_fmt=GENCODE_ANNOT_FMT,
    output:
        GENCODE_GENOME_ANNOT_FILE,
    resources:
        gencode_download_jobs=1,
    retries: config["download"]["retries"]
    wrapper:
        GENCODE_GENOME_ANNOT_WRAPPER
