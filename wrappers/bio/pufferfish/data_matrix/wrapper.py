__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import pandas as pd
import vaex as vx

sample_names = snakemake.params.get("samples")
assert sample_names is not None, "params: samples is a required parameter"

collapse_techreps = snakemake.params.get("collapse_techreps", False)

data_matrix_df = None
for data_file, sample_name in zip(snakemake.input, sample_names):
    if data_file.endswith((".hdf", ".hdf5")):
        data_df = vx.open(data_file)
        data_df.drop(
            [
                data_df.column_names(i)
                for i in range(data_df.shape[1])
                if i not in [0, int(snakemake.params.data_col)]
            ],
            inplace=True,
        )
        data_df.rename(data_df.column_names[-1], sample_name)
        if data_matrix_df:
            if sample_name in data_matrix_df.column_names:
                print(
                    f"Collapsing {snakemake.output[0]} technical "
                    f"replicate {sample_name}",
                    flush=True,
                )
                data_matrix_df[sample_name] = (
                    data_matrix_df[sample_name] + data_df[sample_name]
                )
            else:
                data_matrix_df.join(
                    data_df,
                    how="left",
                    left_on=data_matrix_df.column_names[0],
                    right_on=data_df.column_names[0],
                    allow_duplication=False,
                    cardinality_other=data_df.shape[0],
                    inplace=True,
                )
        else:
            data_matrix_df[sample_name] = data_df
    else:
        data_df = pd.read_csv(
            data_file,
            sep="\t",
            header=None,
            index_col=0,
            usecols=[0, snakemake.params.data_col],
            chunk_size=int(1e7),
            engine="c",
            low_memory=False,
        )
        data_df.columns = [sample_name]
        if data_matrix_df:
            data_matrix_df = pd.concat(
                [data_matrix_df, data_df], axis=1, verify_integrity=True
            )
        else:
            data_matrix_df = pd.DataFrame()
    assert (
        data_matrix_df.shape[0] == data_df.shape[0]
    ), "Quant files do not have same rows"

if isinstance(data_matrix_df, pd.DataFrame):
    data_matrix_df.index.name = "ID_REF"
    if collapse_techreps and data_matrix_df.columns.duplicated().any():
        print(f"Collapsing {snakemake.output[0]} technical replicates", flush=True)
        data_matrix_df = data_matrix_df.groupby(data_matrix_df.columns, axis=1).sum()
    data_matrix_df.sort_index(inplace=True)
    data_matrix_df.to_csv(snakemake.output[0], sep="\t")
else:
    data_matrix_df.rename(data_matrix_df.column_names[0], "ID_REF")
    data_matrix_df = data_matrix_df.sort(by=data_matrix_df.column_names[0])
    data_matrix_df.export_hdf5(
        snakemake.output[0],
        column_count=data_matrix_df.shape[1],
        writer_threads=data_matrix_df.shape[1],
    )
