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


rule fix_gencode_genome_seq_ids:
    input:
        GENCODE_GENOME_SEQ_FILE,
    params:
        pattern=config["gencode"]["seqkit"]["replace"]["search_regex"],
        replacement=lambda w: (
            f"{w.gc_species.title()}_{w.gc_release}_{w.gc_build}_${{1}}"
        ),
    output:
        GENCODE_GENOME_FIXED_SEQ_FILE,
    log:
        GENCODE_GENOME_FIXED_SEQ_LOG,
    threads: config["seqkit"]["replace"]["threads"]
    wrapper:
        SEQKIT_REPLACE_WRAPPER


rule create_gencode_genome_fixed_fasta_list:
    conda:
        "../envs/pandas.yaml"
    input:
        expand(GENCODE_GENOME_FIXED_SEQ_FILE, zip, **EXPAND_PARAMS),
    output:
        GENCODE_GENOME_FIXED_FASTA_LIST_FILE,
    log:
        GENCODE_GENOME_FIXED_FASTA_LIST_LOG,
    script:
        "../scripts/create_file_list_from_paths.py"


rule get_gencode_genome_fixed_seq_ids:
    input:
        GENCODE_GENOME_FIXED_SEQ_FILE,
    params:
        extra="--name",
    output:
        GENCODE_GENOME_FIXED_SEQ_ID_LIST_FILE,
    log:
        GENCODE_GENOME_FIXED_SEQ_ID_LIST_LOG,
    threads: config["seqkit"]["seq"]["threads"]
    wrapper:
        SEQKIT_SEQ_WRAPPER


rule merge_gencode_genome_fixed_seq_ids:
    input:
        expand(GENCODE_GENOME_FIXED_SEQ_ID_LIST_FILE, zip, **EXPAND_PARAMS),
    output:
        GENCODE_GENOME_MERGED_FIXED_SEQ_ID_FILE,
    log:
        GENCODE_GENOME_MERGED_FIXED_SEQ_ID_LOG,
    shell:
        "cat {input} 1> {output} 2> {log}"


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
