__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "BSD 3-Clause"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

infiles = snakemake.input.get("list_file")
if infiles is not None:
    infiles = "--infile-list {infiles}"
else:
    infiles = snakemake.input

pattern = snakemake.params.get("pattern")
extra = snakemake.params.get("extra")
if isinstance(pattern, (list, tuple)):
    flags = []
    if pattern is not None:
        if extra is not None:
            if isinstance(extra, (list, tuple)):
                assert len(extra) == len(
                    pattern
                ), "params: extra must be same length as pattern when list-like"
                for p, e in zip(pattern, extra):
                    flags.append(f"--pattern '{p}' {e}")
        else:
            for p in pattern:
                flags.append(f"--pattern '{p}'")
    elif extra is not None:
        for e in extra:
            flags.append(e)

    cmd_prefix = "seqkit grep --threads {snakemake.threads}"
    shell_cmd = f"{cmd_prefix} {flags[0]} {{infiles}}"
    if len(flags) > 2:
        shell_cmd += " | " + " | ".join(
            f"{cmd_prefix} {f}" for f in flags[1 : len(flags) - 1]
        )
    shell_cmd += (
        f" | {cmd_prefix} {flags[-1]} --out-file {{snakemake.output[0]}} {{log}}"
    )
else:
    flags = ""
    if pattern is not None:
        flags += "--pattern {pattern}"
    id_regexp = snakemake.params.get("id_regexp")
    if id_regexp is not None:
        flags += " --id-regexp '{id_regexp}'"
    extra = snakemake.params.get("extra")
    if extra is not None:
        flags += " {extra}"

    shell(
        "seqkit grep"
        " --threads {snakemake.threads}"
        " {flags}"
        " {infiles}"
        " --out-file {snakemake.output[0]}"
        " {log}"
    )

shell(shell_cmd)
