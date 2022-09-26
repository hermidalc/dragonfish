from os.path import basename, join
import gzip

import vaex as vx
from snakemake.utils import makedirs

makedirs(snakemake.output[0])

file_basename = basename(snakemake.input[0]).partition(".")[0]

split_num = 1
num_split_entries = 0
num_total_entries = 0
split_uniprot_ids = []
split_genbank_ids = []
with gzip.open(snakemake.input[0], "rt") as fh:
    for line in fh:
        uniprot_id, db_name, db_id = line.strip().split("\t")
        if db_name != "EMBL-CDS" or db_id == "-":
            continue
        split_uniprot_ids.append(uniprot_id)
        split_genbank_ids.append(db_id)
        num_split_entries += 1
        num_total_entries += 1
        if num_split_entries == int(float(snakemake.params.split_size)):
            hdf5_file = join(
                snakemake.output[0], f"{file_basename}_{split_num:02}.hdf5"
            )
            print(f"Writing {hdf5_file}")
            vx.from_dict(
                {"genbank_cds_id": split_genbank_ids, "uniprot_id": split_uniprot_ids}
            ).export_hdf5(hdf5_file)
            split_num += 1
            num_split_entries = 0
            split_uniprot_ids = []
            split_genbank_ids = []
print(f"Writing {hdf5_file}")
vx.from_dict(
    {"genbank_cds_id": split_genbank_ids, "uniprot_id": split_uniprot_ids}
).export_hdf5(hdf5_file)
print(f"Parsed {num_total_entries} entries")
