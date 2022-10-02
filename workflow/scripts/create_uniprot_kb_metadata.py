import gzip
from collections import defaultdict
from lxml import etree

import pandas as pd
from Bio.SeqIO.UniprotIO import Parser

num_elems = 0
main_data = defaultdict(list)
dbxref_data = defaultdict(list)
with gzip.open(snakemake.input.kb_file, "rb") as xml_fh:
    for _, elem in etree.iterparse(
        xml_fh, events=("end",), tag="{http://uniprot.org/uniprot}entry"
    ):
        rec = Parser(elem).parse()
        elem.clear()
        del elem
        rec_dbxrefs = defaultdict(list)
        for k, v in [x.split(":", maxsplit=1) for x in rec.dbxrefs]:
            rec_dbxrefs[k.strip()].append(v.strip())
        main_data["id"].append(rec.id)
        main_data["name"].append(rec.name)
        main_data["description"].append(rec.description)
        for db_name in snakemake.params.dbxref_names:
            for db_id in rec_dbxrefs[db_name]:
                dbxref_data["id"].append(rec.id)
                dbxref_data["db"].append(db_name)
                dbxref_data["db_id"].append(db_id)
        num_elems += 1

pd.DataFrame(main_data).to_csv(snakemake.output.main, sep="\t", index=False)
pd.DataFrame(dbxref_data).to_csv(snakemake.output.dbxref, sep="\t", index=False)
