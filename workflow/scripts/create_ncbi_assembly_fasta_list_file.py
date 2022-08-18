from os.path import join, split

import pandas as pd

file_parts = [split(f) for f in snakemake.input]
file_parts_df = pd.DataFrame.from_records(file_parts, columns=["dir", "name"])
file_parts_df.sort_values(
    by=["dir", "name"], ascending=[True, False], inplace=True, ignore_index=True
)
file_parts = file_parts_df.to_records(index=False).tolist()
files = [join(*p) for p in file_parts]
with open(snakemake.output[0], "w") as fh:
    fh.write("{}\n".format("\n".join(files)))
