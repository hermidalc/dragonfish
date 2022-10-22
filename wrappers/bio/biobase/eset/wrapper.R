# __author__ = "Leandro C. Hermida"
# __email__ = "hermidalc@pitt.edu"
# __license__ = "BSD 3-Clause"

suppressPackageStartupMessages({
    library(Biobase)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

adata <- read.delim(
    snakemake@input[["assay"]], sep = "\t", header = TRUE, row.names = 1
)
pdata <- read.delim(
    snakemake@input[["pheno"]], sep = "\t", header = TRUE,
    row.names = snakemake@params[["samples"]]
)
fdata <- read.delim(
    snakemake@input[["annot"]], sep = "\t", header = TRUE, row.names = 1
)

eset <- ExpressionSet(
    assayData = as.matrix(adata),
    phenoData = AnnotatedDataFrame(pdata),
    featureData = AnnotatedDataFrame(fdata),
)

saveRDS(eset, file = snakemake@output[[1]])

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()