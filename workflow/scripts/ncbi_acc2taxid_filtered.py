import gzip

import pandas as pd

summary_df = pd.read_csv(
    snakemake.input[0],
    sep="\t",
    index_col="assembly_accession",
    engine="c",
    low_memory=False,
)
with open(snakemake.output[0], "wt") as out_fh:
    for i, gz_file in enumerate(sorted(snakemake.input.files, reverse=True)):
        with gzip.open(gz_file, "rt") as in_fh:
            header = in_fh.readline()
            if i == 0:
                out_fh.write(header)
            for line in in_fh:
                tax_id = line.split("\t")[2].strip()
                if tax_id in summary_df["taxid"]:
                    out_fh.write(line)
