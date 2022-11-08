import pandas as pd

summary_df = pd.read_csv(
    snakemake.input[0],
    sep="\t",
    index_col="assembly_accession",
    engine="c",
    low_memory=False,
)
with open(snakemake.output[0], "wt") as out_fh:
    with open(snakemake.input.nodes, "rt") as in_fh:
        for line in in_fh:
            fields = line.split("|")
            fields = [f.strip() for f in fields]
            if fields[0] in summary_df["taxid"]:
                out_fh.write(f"{'\t|\t'.join(fields)}\t|\n")
