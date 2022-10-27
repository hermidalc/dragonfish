import gzip

with open(snakemake.output[0], "wt") as out_fh:
    for i, gz_file in enumerate(sorted(snakemake.input, reverse=True)):
        with gzip.open(gz_file, "rt") as in_fh:
            if i == 0:
                out_fh.write(in_fh.readline())
            out_fh.write(in_fh.readline())
