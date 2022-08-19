rule get_gencode_genome_seq:
    params:
        protocol=GENCODE_PROTOCOL,
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
        regions=GENCODE_REGIONS,
    output:
        GENCODE_GENOME_SEQ_FILE,
    wrapper:
        GENCODE_GENOME_SEQ_WRAPPER


rule get_gencode_genome_annot:
    params:
        protocol=GENCODE_PROTOCOL,
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        regions=GENCODE_REGIONS,
        annot_fmt=GENCODE_ANNOT_FMT,
    output:
        GENCODE_GENOME_ANNOT_FILE,
    wrapper:
        GENCODE_GENOME_ANNOT_WRAPPER


rule fix_gencode_genome_seq_ids:
    input:
        GENCODE_GENOME_SEQ_FILE,
    params:
        pattern=GENCODE_SEQKIT_REPLACE_SEARCH_REGEX,
        replacement=lambda wildcards: f"${{1}}_{wildcards.gc_build}",
        extra=config["seqkit"]["replace"]["extra_params"],
    output:
        GENCODE_GENOME_FIXED_SEQ_FILE,
    log:
        GENCODE_GENOME_FIXED_SEQ_LOG,
    threads: config["seqkit"]["replace"]["threads"]
    wrapper:
        SEQKIT_REPLACE_WRAPPER


rule merge_gencode_genome_seq:
    input:
        expand(GENCODE_GENOME_FIXED_SEQ_FILE, zip, **EXPAND_PARAMS),
    output:
        GENCODE_GENOME_MERGED_SEQ_FILE,
    log:
        GENCODE_GENOME_MERGED_SEQ_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"
