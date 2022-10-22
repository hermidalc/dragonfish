__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

import pandas as pd

sample_names = snakemake.params.get("samples")
assert sample_names is not None, "params: samples is a required parameter"

data_matrix_df = pd.DataFrame()
for quant_file, sample_name, strand in zip(snakemake.input[0], sample_names):
    quants = pd.read_csv(quant_file, sep="\t", header=None, index_col=0, usecols=[0])
    quants.columns = [sample_name]
    data_matrix_df = pd.concat([data_matrix_df, quants], axis=1, verify_integrity=True)
    assert (
        data_matrix_df.shape[0] == quants.shape[0]
    ), "Quant files do not have same rows"

data_matrix_df.index.name = "ID_REF"

collapse_techreps = snakemake.params.get("collapse_techreps", False)
if collapse_techreps and data_matrix_df.columns.duplicated().any():
    print(f"Collapsing {snakemake.output[0]} technical replicates", flush=True)
    data_matrix_df = data_matrix_df.groupby(data_matrix_df.columns, axis=1).sum()

data_matrix_df.sort_index(inplace=True)
data_matrix_df.to_csv(snakemake.output[0], sep="\t")
