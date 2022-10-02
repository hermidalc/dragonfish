import gzip
from collections import defaultdict
from operator import itemgetter

import pandas as pd
from Bio.SeqIO.UniprotIO import Parser
from lxml import etree

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
        main_data["id"].append(rec.id)
        main_data["name"].append(rec.name)
        main_data["description"].append(rec.description)
        for db_name, db_id in sorted(
            [x.split(":", maxsplit=1) for x in rec.dbxrefs], key=itemgetter(0)
        ):
            dbxref_data["id"].append(rec.id)
            dbxref_data["db"].append(db_name.strip())
            dbxref_data["db_id"].append(db_id.strip())
        num_elems += 1

pd.DataFrame(main_data).to_csv(snakemake.output.main, sep="\t", index=False)
pd.DataFrame(dbxref_data).to_csv(snakemake.output.dbxref, sep="\t", index=False)
