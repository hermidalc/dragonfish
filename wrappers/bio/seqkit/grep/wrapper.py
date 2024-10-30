__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

infiles = snakemake.input.get("list_file")
if infiles:
    infiles = f"--infile-list {infiles}"
else:
    infiles = snakemake.input

pattern = snakemake.params.get("pattern")
extra = snakemake.params.get("extra")
if isinstance(pattern, (list, tuple)):
    flags = []
    if pattern:
        if extra:
            if isinstance(extra, (list, tuple)):
                assert len(extra) == len(
                    pattern
                ), "params: extra must be same length as pattern when list-like"
                for p, e in zip(pattern, extra):
                    flags.append(f"--pattern '{p}' {e}")
        else:
            for p in pattern:
                flags.append(f"--pattern '{p}'")
    elif extra:
        for e in extra:
            flags.append(e)

    cmd_prefix = "seqkit grep --threads {snakemake.threads}"
    shell_cmd = f"{cmd_prefix} {flags[0]} {{infiles}}"
    if len(flags) > 1:
        shell_cmd += " | " + " | ".join(f"{cmd_prefix} {f}" for f in flags[1:])
    shell_cmd += " --out-file {snakemake.output[0]} {log}"
    shell(shell_cmd)
else:
    flags = ""
    if pattern:
        flags += f"--pattern '{pattern}'"
    id_regexp = snakemake.params.get("id_regexp")
    if extra:
        flags += f" {extra}"
    if id_regexp:
        flags += f" --id-regexp '{id_regexp}'"

    shell(
        "seqkit grep"
        " --threads {snakemake.threads}"
        " {flags}"
        " {infiles}"
        " --out-file {snakemake.output[0]}"
        " {log}"
    )
