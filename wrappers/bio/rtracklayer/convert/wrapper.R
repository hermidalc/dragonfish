# __author__ = "Leandro C. Hermida"
# __email__ = "hermidalc@pitt.edu"
# __license__ = "BSD 3-Clause"

suppressPackageStartupMessages({
    library(rtracklayer)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

annot <- import(snakemake@output[[1]])
export(annot, snakemake@output[[1]], format = snakemake@params[["format"]])

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
