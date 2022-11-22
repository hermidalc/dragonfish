import gzip

import pandas as pd

summary_df = pd.read_csv(
    snakemake.input.summary,
    sep="\t",
    index_col="assembly_accession",
    engine="c",
    low_memory=False,
)

tax_ids = set()
summary_tax_ids = set(summary_df["taxid"].values)
with open(snakemake.output[0], "wt") as out_fh:
    for i, gz_file in enumerate(sorted(snakemake.input.files, reverse=True)):
        with gzip.open(gz_file, "rt") as in_fh:
            header = in_fh.readline()
            for line in in_fh:
                fields = line.split("\t")
                tax_id = int(fields[2].strip())
                if tax_id in summary_tax_ids:
                    tax_ids.add(tax_id)
                    ref_id = fields[1].strip()
                    out_fh.write(f"{ref_id}\t{tax_id}\n")

missing_taxids = list(sorted(summary_tax_ids.difference(tax_ids)))
if missing_taxids:
    print(f"Not all assembly taxids were found in acc2taxid: {missing_taxids}")
