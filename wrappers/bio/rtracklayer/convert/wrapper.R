# __author__ = "Leandro C. Hermida"
# __email__ = "hermidalc@pitt.edu"
# __license__ = "MIT"

suppressPackageStartupMessages({
    library(rtracklayer)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

gzip <- endsWith(tolower(snakemake@output[[1]]), ".gz")
if (gzip) {
    export_file <- tools::file_path_sans_ext(snakemake@output[[1]])
} else {
    export_file <- snakemake@output[[1]]
}

export(
    import(snakemake@input[[1]]),
    con = export_file,
    format = snakemake@params[["format"]]
)

if (gzip) {
    R.utils::gzip(
        export_file,
        destname = snakemake@output[[1]],
        overwrite = TRUE, skip = FALSE, remove = TRUE
    )
}

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
