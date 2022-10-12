__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

cmd = snakemake.params.get("cmd")
assert cmd is not None, "input: cmd is a required input parameter"

infiles = snakemake.input.get("list_file")
if infiles is not None:
    infiles = f"--infile-list {infiles}"
else:
    infiles = snakemake.input

pattern = snakemake.params.get("pattern")
extra = snakemake.params.get("extra")
if cmd == "grep" and isinstance(pattern, (list, tuple)):
    flags = []
    if pattern is not None:
        if extra is not None:
            if isinstance(extra, (list, tuple)):
                assert len(extra) == len(
                    pattern
                ), "extra must be same length as pattern when list-like"
                for p, e in zip(pattern, extra):
                    flags.append([f"--pattern '{p}' {e}"])
        else:
            for p in pattern:
                flags.append([f"--pattern '{p}'"])
    elif extra is not None:
        for e in extra:
            flags.append([e])

    cmd_prefix = f"seqkit {cmd} --threads {snakemake.threads}"
    shell_cmd = f" {flags[0]} {infiles} |"
    shell_cmd += " | ".join(f"{cmd_prefix} {f}" for f in flags[1 : len(flags)])
    shell_cmd += f"{cmd_prefix} --out-file {snakemake.output[0]} {log}"
else:
    flags = ""
    if pattern is not None:
        flags += f"--pattern '{pattern}'"
        replacement = snakemake.params.get("replacement")
        if replacement is not None:
            flags += f" --replacement '{replacement}'"
    if extra is not None:
        flags += f" {extra}"

    shell_cmd = (
        f"seqkit {cmd}"
        f" --threads {snakemake.threads}"
        f" {flags}"
        f" {infiles}"
        f" --out-file {snakemake.output[0]}"
        f" {log}"
    )

shell(shell_cmd)
