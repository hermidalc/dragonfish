import gzip
import re
from lxml import etree
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


def write_split_file(elems, out_dir, file_basename, num, sep):
    split_filename = f"{file_basename}_{num:04}.xml.gz"
    file = join(out_dir, split_filename)
    makedirs(out_dir)
    print(f"Writing {split_filename}")
    with gzip.open(file, "wt") as fh:
        fh.write(xml_header)
        fh.write(sep.join(elems))
        fh.write(sep + xml_footer)


split_num = 1
split_elems = []
num_split_elems = 0
num_total_elems = 0
if snakemake.params.parser == "python_lxml":
    elem_sep = ""
    open_mode = "rb"
else:
    elem_sep = "\n"
    open_mode = "rt"
with gzip.open(snakemake.input[0], open_mode) as xml_fh:
    if snakemake.params.parser == "python_lxml":
        for _, elem in etree.iterparse(
            xml_fh, events=("end",), tag=f"{{{xmlns}}}entry"
        ):
            split_elems.append(etree.tostring(elem, encoding="unicode"))
            elem.clear()
            del elem
            num_split_elems += 1
            num_total_elems += 1
            if num_split_elems == int(float(snakemake.params.split_size)):
                write_split_file(
                    split_elems,
                    snakemake.output[0],
                    snakemake.params.basename,
                    split_num,
                    elem_sep,
                )
                split_num += 1
                split_elems = []
                num_split_elems = 0
    else:
        in_elem = False
        elem_lines = []
        elem_s_regex = re.compile(r"\s*<\s*entry\s+.*?>")
        elem_e_regex = re.compile(r"\s*</\s*entry\s*>")
        for line in xml_fh:
            line = line.rstrip()
            if elem_s_regex.match(line):
                in_elem = True
                elem_lines.append(line)
            elif elem_e_regex.match(line):
                assert (
                    in_elem
                ), "Closing entry tag without previously reading opening tag"
                elem_lines.append(line)
                in_elem = False
                split_elems.append("\n".join(elem_lines))
                elem_lines = []
                num_split_elems += 1
                num_total_elems += 1
                if num_split_elems == int(float(snakemake.params.split_size)):
                    write_split_file(
                        split_elems,
                        snakemake.output[0],
                        snakemake.params.basename,
                        split_num,
                        elem_sep,
                    )
                    split_num += 1
                    split_elems = []
                    num_split_elems = 0
            elif in_elem:
                elem_lines.append(line)
write_split_file(
    split_elems, snakemake.output[0], snakemake.params.basename, split_num, elem_sep
)
print(f"Parsed {num_total_elems} {snakemake.params.basename} records")
