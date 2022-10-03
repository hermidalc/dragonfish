import gzip
import re
from collections import defaultdict
from operator import itemgetter

import pandas as pd
from Bio.SeqIO.UniprotIO import Parser
from lxml import etree

proteome_df = pd.read_csv(
    snakemake.input.proteome_file, sep="\t", index_col="Proteome Id"
)

split_offset = list(
    range(0, snakemake.params.kb_size, int(float(snakemake.params.split_size)))
)[int(snakemake.params.split_num) - 1]

entry_lines = []
entry_offset = 0
parse_entries = False
num_entries_parsed = 0
entry_s_regex = re.compile(r"^\s*<\s*entry\s+.*?>")
entry_e_regex = re.compile(r"^\s*</\s*entry\s*>")
dbxrefs = defaultdict(list)
with gzip.open(snakemake.input.kb_file, "rt") as xml_fh:
    for line in xml_fh:
        if entry_s_regex.match(line):
            if entry_offset == split_offset:
                parse_entries = True
            entry_offset += 1
        if parse_entries:
            entry_lines.append(line)
            if entry_e_regex.match(line):
                entry = "".join(entry_lines)
                elem = etree.fromstring(entry)
                rec = Parser(elem).parse()
                entry_lines = []
                rec_dbxrefs = defaultdict(list)
                for k, v in sorted(
                    [x.split(":", maxsplit=1) for x in rec.dbxrefs], key=itemgetter(0)
                ):
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
                num_entries_parsed += 1
                if num_entries_parsed == int(float(snakemake.params.split_size)):
                    break

pd.DataFrame(dbxrefs).to_csv(snakemake.output[0], sep="\t", index=False)
