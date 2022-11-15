import gzip

from gffutils import create_db, Feature

db = create_db(
    snakemake.input[0], dbfn=":memory:", keep_order=True, merge_strategy="create_unique"
)

with gzip.open(snakemake.output[0], "wt") if snakemake.output[0].endswith(
    ".gz"
) else open(snakemake.output[0], "wt") as out_fh:
    for d in db.directives:
        out_fh.write(f"##{d}\n")
    for parent in db.features_of_type(("gene", "mRNA")):
        if parent.featuretype == "gene" and list(
            db.children(parent, featuretype="mRNA")
        ):
            continue
        if "protein_id" in parent.attributes:
            out_fh.write(f"{parent}\n")
            for exon in db.children(parent, featuretype="exon"):
                exon.attributes["protein_id"] = parent.attributes["protein_id"]
                out_fh.write(f"{exon}\n")
        else:
            cdses, protein_ids = [], []
            for c in db.children(parent, featuretype="CDS", order_by="start"):
                for a in c.attributes:
                    if a == "protein_id":
                        cdses.append(c)
                        protein_ids.extend(c.attributes["protein_id"])
            protein_ids = list(set(protein_ids))
            if len(protein_ids) > 1:
                # XXX: hack to fix edge case
                # e.g. GCA_002221805.1_ViralProj14310 ID=gene-cI
                print(
                    f"Warning: {snakemake.input[0]}: {parent.id}:"
                    " gene with no mRNA(s) multi protein_id fix"
                )
                exons = [
                    e for e in db.children(parent, featuretype="exon", order_by="start")
                ]
                for cds in cdses:
                    no_cds_exon = True
                    for exon in exons:
                        if (
                            no_cds_exon
                            and "protein_id" not in exon.attributes
                            and cds.start >= exon.start
                            and cds.end <= exon.end
                        ):
                            exon.attributes["protein_id"] = cds.attributes["protein_id"]
                            no_cds_exon = False
                    if no_cds_exon:
                        exons.append(
                            Feature(
                                seqid=cds.seqid,
                                source=cds.source,
                                featuretype="exon",
                                start=cds.start,
                                end=cds.end,
                                score=".",
                                strand=cds.strand,
                                frame=".",
                                attributes={
                                    "protein_id": cds.attributes["protein_id"],
                                    "Parent": [parent.id],
                                },
                            )
                        )
                out_fh.write(f"{parent}\n")
                for exon in exons:
                    out_fh.write(f"{exon}\n")
            elif len(protein_ids) == 1:
                parent.attributes["protein_id"] = protein_ids
                out_fh.write(f"{parent}\n")
                for exon in db.children(parent, featuretype="exon"):
                    exon.attributes["protein_id"] = parent.attributes["protein_id"]
                    out_fh.write(f"{exon}\n")
            else:
                print(
                    f"Warning: {snakemake.input[0]}: {parent.id}"
                    " no child CDS protein_id attributes found"
                )
