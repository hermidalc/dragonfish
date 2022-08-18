from os.path import join, split

import pandas as pd

parts = [split(f) for f in snakemake.input]
parts_df = pd.DataFrame.from_records(parts, columns=["dir", "name"])
parts_df.sort_values(
    by=["dir", "name"],
    ascending=[True, False],
    inplace=True,
    ignore_index=True,
)
parts = parts_df.to_records(index=False).tolist()
files = [join(*p) for p in parts]
with open(snakemake.output[0], "w") as fh:
    fh.write("{}\n".format("\n".join(files)))
