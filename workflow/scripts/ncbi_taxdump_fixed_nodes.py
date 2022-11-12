# XXX: mappings for missing valid ranks in Cedar
rank_map = {
    "biotype": "no rank",
    "clade": "class",
    "forma specialis": "subspecies",
    "genotype": "no rank",
    "isolate": "no rank",
    "morph": "forma",
    "pathogroup": "no rank",
    "section": "subgenus",
    "subsection": "subgenus",
    "series": "species group",
    "serotype": "no rank",
    "serogroup": "no rank",
    "strain": "no rank",
    "subcohort": "cohort",
}

with open(snakemake.output[0], "wt") as out_fh:
    with open(snakemake.input[0], "rt") as in_fh:
        for line in in_fh:
            fields = [f.strip() for f in line.split("|")]
            tax_id, parent_tax_id, rank = fields[0:3]
            rank = rank.lower()
            if rank in rank_map:
                rank = rank_map[rank]
            out_line = "\t|\t".join([tax_id, parent_tax_id, rank] + fields[3:])
            out_fh.write(f"{out_line}\t|\n")
