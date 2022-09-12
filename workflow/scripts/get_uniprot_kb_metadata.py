import gzip
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from joblib import delayed, Parallel


def parse_xml_file(file, proteome_df, dbxref_names):
    file_main_data = defaultdict(list)
    file_dbxref_data = defaultdict(list)
    with gzip.open(file, "rt") as xml_fh:
        for rec in SeqIO.parse(xml_fh, "uniprot-xml"):
            rec_dbxrefs = defaultdict(list)
            for k, v in [x.split(":", maxsplit=1) for x in rec.dbxrefs]:
                rec_dbxrefs[k.strip()].append(v.strip())
            if (
                "Proteomes" not in rec_dbxrefs
                or not pd.Series(rec_dbxrefs["Proteomes"]).isin(proteome_df.index).any()
                or not pd.Series(rec_dbxrefs.keys()).isin(dbxref_names).any()
            ):
                continue
            file_main_data["id"].append(rec.id)
            file_main_data["name"].append(rec.name)
            file_main_data["description"].append(rec.description)
            for db_name in dbxref_names:
                for db_id in rec_dbxrefs[db_name]:
                    file_dbxref_data["id"].append(rec.id)
                    file_dbxref_data["db"].append(db_name)
                    file_dbxref_data["db_id"].append(db_id)
    return file_main_data, file_dbxref_data


proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

main_data, dbxref_data = zip(
    *Parallel(
        n_jobs=snakemake.threads,
        backend=snakemake.params.backend,
        verbose=snakemake.params.verbosity,
    )(
        delayed(parse_xml_file)(file, proteome_df, snakemake.params.dbxref_names)
        for file in sorted(snakemake.input.kb_file)
    )
)

main_data = {
    k1: [v for d in main_data for k2, l in d.items() for v in l if k2 == k1]
    for k1 in main_data[0].keys()
}
dbxref_data = {
    k1: [v for d in dbxref_data for k2, l in d.items() for v in l if k2 == k1]
    for k1 in dbxref_data[0].keys()
}

pd.DataFrame(main_data).to_csv(snakemake.output.main, sep="\t", index=False)
pd.DataFrame(dbxref_data).to_csv(snakemake.output.dbxref, sep="\t", index=False)
