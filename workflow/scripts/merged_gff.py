import gzip

with gzip.open(snakemake.output[0], "wt") if snakemake.output[0].endswith(
    ".gz"
) else open(snakemake.output[0], "wt") as out_fh:
    out_fh.write(f"##gff-version 3\n")
    for gtf_file in snakemake.input:
        with gzip.open(gtf_file, "rt") if gtf_file.endswith(".gz") else open(
            gtf_file, "rt"
        ) as in_fh:
            lines_written = False
            for line in in_fh:
                if not line.startswith("#"):
                    out_fh.write(line)
                    lines_written = True
            if lines_written:
                out_fh.write("###\n")
