import gzip
import re
from os.path import join

from snakemake.utils import makedirs

xmlns = "http://uniprot.org/uniprot"

xml_header = f"""<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="{xmlns}"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">
"""

xml_footer = f"""<copyright>
Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms Distributed under the Creative Commons Attribution (CC BY 4.0) License
</copyright>
</uniprot>
"""


def write_split_file(entries, out_dir, file_basename, num):
    split_filename = f"{file_basename}_{num:04}.xml.gz"
    file = join(out_dir, split_filename)
    makedirs(out_dir)
    print(f"Writing {split_filename}")
    with gzip.open(file, "wt") as fh:
        fh.write(xml_header)
        fh.write("".join(entries))
        fh.write("" + xml_footer)


split_num = 1
entry_lines = []
num_split_entries = 0
num_total_entries = 0
in_entry = False
entry_s_regex = re.compile(r"\s*<\s*entry\s+.*?>")
entry_e_regex = re.compile(r"\s*</\s*entry\s*>")
with gzip.open(snakemake.input[0], "rt") as xml_fh:
    for line in xml_fh:
        if entry_s_regex.match(line):
            in_entry = True
            entry_lines.append(line)
        elif entry_e_regex.match(line):
            assert in_entry, "Closing entry tag without previously reading opening tag"
            in_entry = False
            entry_lines.append(line)
            num_split_entries += 1
            num_total_entries += 1
            if num_split_entries == int(float(snakemake.params.split_size)):
                write_split_file(
                    entry_lines,
                    snakemake.output[0],
                    snakemake.params.basename,
                    split_num,
                )
                entry_lines = []
                num_split_entries = 0
                split_num += 1
        elif in_entry:
            entry_lines.append(line)
write_split_file(entry_lines, snakemake.output[0], snakemake.params.basename, split_num)
print(f"Parsed {num_total_entries} {snakemake.params.basename} records")
