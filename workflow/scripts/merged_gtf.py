import gzip

with gzip.open(snakemake.output[0], "wt") as out_fh:
    out_fh.write("##gff-version 2\n")
    for gtf_file in snakemake.input:
        with gzip.open(gtf_file, "rt") as in_fh:
            for line in in_fh:
                if not line.startswith("#"):
                    out_fh.write(line)
            out_fh.write("###\n")
