rule get_gencode_genome_fasta:
    params:
        protocol=GENCODE_PROTOCOL,
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
        regions=GENCODE_REGIONS,
    output:
        GENCODE_GENOME_FASTA_FILE,
    wrapper:
        GENCODE_GENOME_SEQ_WRAPPER


rule fix_gencode_genome_fasta_ids:
    input:
        GENCODE_GENOME_FASTA_FILE,
    params:
        cmd="replace",
        pattern=config["gencode"]["seqkit"]["replace"]["search_regex"],
        replacement=lambda w: f"{w.gc_species.title()}_{w.gc_release}_{w.gc_build}_${{1}}",
        extra=f'--line-width {config["seqkit"]["line_width"]}',
    output:
        GENCODE_GENOME_FIXED_FASTA_FILE,
    log:
        GENCODE_GENOME_FIXED_FASTA_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule get_gencode_genome_fixed_fasta_ids:
    input:
        GENCODE_GENOME_FIXED_FASTA_FILE,
    params:
        cmd="seq",
        extra="--name",
    output:
        GENCODE_GENOME_FIXED_FASTA_ID_LIST_FILE,
    log:
        GENCODE_GENOME_FIXED_FASTA_ID_LIST_LOG,
    threads: config["seqkit"]["threads"]
    wrapper:
        SEQKIT_WRAPPER


rule merge_gencode_genome_fixed_fasta_ids:
    input:
        expand(GENCODE_GENOME_FIXED_FASTA_ID_LIST_FILE, zip, **EXPAND_PARAMS),
    output:
        GENCODE_GENOME_MERGED_FIXED_FASTA_ID_FILE,
    log:
        GENCODE_GENOME_MERGED_FIXED_FASTA_ID_LOG,
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
