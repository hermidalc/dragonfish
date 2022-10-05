import gzip

split_offsets = list(
    range(0, int(snakemake.params.kb_size), int(float(snakemake.params.split_size)))
)

entry_offset = 0
split_pos = []
with gzip.open(snakemake.input[0], "rb") as kb_fh:
    while line := kb_fh.readline():
        if split_offsets and line.decode("utf-8").startswith("<entry "):
            if entry_offset == split_offsets[0]:
                split_pos.append(kb_fh.tell() - len(line))
                split_offsets.pop(0)
            entry_offset += 1

with open(snakemake.output[0], "wt") as out_fh:
    out_fh.writelines(f"{p}\n" for p in split_pos)
