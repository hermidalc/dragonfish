import gzip
from os.path import basename, join
from tempfile import gettempdir, TemporaryDirectory

import vaex as vx

split_num = 1
num_split_entries = 0
num_total_entries = 0
split_uniprot_ids = []
split_genbank_ids = []
file_basename = basename(snakemake.input[0]).partition(".")[0]
with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    with gzip.open(snakemake.input[0], "rt") as fh:
        for line in fh:
            uniprot_id, db_name, db_id = line.strip().split("\t")
            if db_name == "EMBL-CDS" and db_id != "-":
                split_uniprot_ids.append(uniprot_id)
                split_genbank_ids.append(db_id)
                num_split_entries += 1
                num_total_entries += 1
                if num_split_entries == int(float(snakemake.params.split_size)):
                    idmap_file = join(tmp_dir, f"{file_basename}_{split_num:02}.hdf5")
                    print(f"Writing {idmap_file}")
                    vx.from_dict(
                        {
                            "genbank_id": split_genbank_ids,
                            "uniprot_id": split_uniprot_ids,
                        }
                    ).export_hdf5(idmap_file, column_count=2, writer_threads=2)
                    split_num += 1
                    num_split_entries = 0
                    split_uniprot_ids = []
                    split_genbank_ids = []
    idmap_file = join(tmp_dir, f"{file_basename}_{split_num:02}.hdf5")
    print(f"Writing {idmap_file}")
    vx.from_dict(
        {"genbank_id": split_genbank_ids, "uniprot_id": split_uniprot_ids}
    ).export_hdf5(idmap_file, column_count=2, writer_threads=2)
    vx.open(join(tmp_dir, f"{file_basename}_*.hdf5")).export_hdf5(
        snakemake.output[0], column_count=2, writer_threads=2
    )

print(f"Parsed {num_total_entries} entries")
