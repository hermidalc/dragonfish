import gzip

with (
    gzip.open(snakemake.output[0], "wt")
    if snakemake.output[0].endswith(".gz")
    else open(snakemake.output[0], "wt")
) as out_fh:
    out_fh.write(snakemake.params.header)
    for gtf_file in snakemake.input:
        with (
            gzip.open(gtf_file, "rt")
            if gtf_file.endswith(".gz")
            else open(gtf_file, "rt")
        ) as in_fh:
            for line in in_fh:
                if not line.startswith("#"):
                    out_fh.write(line)
