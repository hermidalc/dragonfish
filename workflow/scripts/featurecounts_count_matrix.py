from os.path import basename, join, splitext
from tempfile import gettempdir, TemporaryDirectory

import pandas as pd
import vaex as vx

with TemporaryDirectory(dir=snakemake.resources.get("tmpdir", gettempdir())) as tmp_dir:
    file_basename = basename(snakemake.input[0]).partition(".")[0]
    cols = pd.read_csv(snakemake.input[0], sep="\t", nrows=1, comment="#").columns
    use_cols = [0] + list(range(6, len(cols)))
    col_names = ["genbank_id"] + [splitext(basename(c))[0] for c in cols[use_cols[1:]]]
    for i, df in enumerate(
        vx.read_csv(
            snakemake.input[0],
            sep="\t",
            chunk_size=int(1e7),
            engine="c",
            low_memory=False,
            comment="#",
            header=0,
            usecols=use_cols,
            names=col_names,
        ),
        start=1,
    ):
        df.export_hdf5(join(tmp_dir, f"{file_basename}_{i:02}.hdf5"))
    df = vx.open(join(tmp_dir, f"{file_basename}_*.hdf5"))
    df.export_hdf5(
        snakemake.output[0], column_count=df.shape[1], writer_threads=df.shape[1]
    )
