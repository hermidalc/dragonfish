import gzip

with open(snakemake.output[0], "wt") as out_fh:
    for i, gz_file in enumerate(sorted(snakemake.input.files, reverse=True)):
        with gzip.open(gz_file, "rt") as in_fh:
            if i == 0:
                out_fh.write(in_fh.readline())
            line = in_fh.readline()
            tax_id = line.split("\t")[2].strip()
            out_fh.write(line)
