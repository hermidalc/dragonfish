import gzip
from collections import defaultdict

import pandas as pd
from Bio.SeqIO.UniprotIO import Parser as UniprotParser
from lxml import etree

proteome_df = pd.read_csv(snakemake.input.proteome, sep="\t", index_col="Proteome Id")

with open(snakemake.input.split_pos, "rt") as pos_fh:
    split_pos = pos_fh.read().splitlines()

split_num = int(snakemake.params.split_num)
start_pos = int(split_pos[split_num - 1])
end_pos = int(split_pos[split_num]) if split_num < len(split_pos) else -1

entry_lines = []
dbxrefs = defaultdict(list)
with gzip.open(snakemake.input.kb, "rb") as kb_fh:
    kb_fh.seek(start_pos)
    while line := kb_fh.readline():
        entry_lines.append(line)
        if line.decode("utf-8").startswith("</entry>"):
            rec = UniprotParser(etree.fromstring(b"".join(entry_lines))).parse()
            entry_lines = []
            rec_dbxrefs = defaultdict(list)
            for k, v in [x.split(":", maxsplit=1) for x in sorted(rec.dbxrefs)]:
                rec_dbxrefs[k.strip()].append(v.strip())
            if (
                "Proteomes" in rec_dbxrefs
                and any(x in proteome_df.index for x in rec_dbxrefs["Proteomes"])
                and any(x in rec_dbxrefs for x in snakemake.params.dbxref_names)
            ):
                for db_name in snakemake.params.dbxref_names:
                    for db_id in rec_dbxrefs[db_name]:
                        dbxrefs["id"].append(rec.id)
                        dbxrefs["db"].append(db_name)
                        dbxrefs["db_id"].append(db_id)
            if kb_fh.tell() == end_pos:
                break

pd.DataFrame(dbxrefs).to_csv(snakemake.output[0], sep="\t", index=False)
