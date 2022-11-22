import gzip

with open(snakemake.output[0], "wt") as out_fh:
    with gzip.open(snakemake.input[0], "rt") as in_fh:
        for line in in_fh:
            uniprot_id, db_name, db_id = line.strip().split("\t")
            if db_name == "NCBI_TaxID" and db_id != "-":
                # Write zero (0) for last gi column
                out_fh.write(f"{uniprot_id}\t{uniprot_id}\t{db_id}\t0\n")
