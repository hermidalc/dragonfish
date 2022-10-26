from os.path import basename, join
from tempfile import gettempdir, TemporaryDirectory

import vaex as vx

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    file_basename = basename(snakemake.input[0]).partition(".")[0]
    for i, df in enumerate(
        vx.read_csv(snakemake.input[0], chunk_size=snakemake.params.chunk_size), start=1
    ):
        df.export_hdf5(join(tmp_dir, f"{file_basename}_{i:03}.hdf5"))
    vx.open(join(tmp_dir, f"{file_basename}_*.hdf5")).export_hdf5(
        snakemake.output[0], column_count=3, writer_threads=3
    )
