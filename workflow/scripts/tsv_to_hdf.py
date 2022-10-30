from os.path import basename, join
from tempfile import gettempdir, TemporaryDirectory

import vaex as vx

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    file_basename = basename(snakemake.input[0]).partition(".")[0]
    for i, df in enumerate(
        vx.read_csv(
            snakemake.input[0],
            sep="\t",
            chunk_size=int(1e7),
            engine="c",
            low_memory=False,
        ),
        start=1,
    ):
        df.export_hdf5(join(tmp_dir, f"{file_basename}_{i:03}.hdf5"))
    df = vx.open(join(tmp_dir, f"{file_basename}_*.hdf5"))
    df.export_hdf5(
        snakemake.output[0], column_count=df.shape[1], writer_threads=df.shape[1]
    )
