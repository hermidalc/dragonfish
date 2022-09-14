import gzip
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

main_data = defaultdict(list)
dbxref_data = defaultdict(list)
with gzip.open(snakemake.input.kb_file, "rt") as xml_fh:
    for rec in SeqIO.parse(xml_fh, "uniprot-xml"):
        rec_dbxrefs = defaultdict(list)
        for k, v in [x.split(":", maxsplit=1) for x in rec.dbxrefs]:
            rec_dbxrefs[k.strip()].append(v.strip())
        if (
            "Proteomes" not in rec_dbxrefs
            or not pd.Series(rec_dbxrefs["Proteomes"]).isin(proteome_df.index).any()
            or not pd.Series(rec_dbxrefs.keys())
            .isin(snakemake.params.dbxref_names)
            .any()
        ):
            continue
        main_data["id"].append(rec.id)
        main_data["name"].append(rec.name)
        main_data["description"].append(rec.description)
        for db_name in snakemake.params.dbxref_names:
            for db_id in rec_dbxrefs[db_name]:
                dbxref_data["id"].append(rec.id)
                dbxref_data["db"].append(db_name)
                dbxref_data["db_id"].append(db_id)

pd.DataFrame(main_data).to_csv(snakemake.output.main, sep="\t", index=False)
pd.DataFrame(dbxref_data).to_csv(snakemake.output.dbxref, sep="\t", index=False)
