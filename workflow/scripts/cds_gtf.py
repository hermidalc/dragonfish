import gzip

import gffutils

db = gffutils.create_db(
    snakemake.input[0], dbfn=":memory:", keep_order=True, merge_strategy="create_unique"
)

with gzip.open(snakemake.output[0], "wt") if snakemake.output[0].endswith(
    ".gz"
) else open(snakemake.output[0], "wt") as out_fh:
    out_fh.write(f"#gtf-version 2.2\n")
    for cds in db.features_of_type("CDS"):
        if "protein_id" in cds.attributes:
            assert len(cds.attributes["protein_id"]) == 1, (
                f"Warning: {snakemake.input[0]}: {cds}: more than one protein_id"
                " attribute"
            )
            protein_id = cds.attributes["protein_id"][0]
            cds_parts = str(cds).split("\t")
            cds_parts[8] = f'protein_id "{protein_id}";'
            cds_str = "\t".join(cds_parts)
            out_fh.write(f"{cds_str}\n")
        elif "pseudo" not in cds.attributes or any(
            p != "true" for p in cds.attributes["pseudo"]
        ):
            print(f"Warning: {snakemake.input[0]}: {cds}: no protein_id attribute")
