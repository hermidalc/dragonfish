__author__ = "Leandro C. Hermida"
__email__ = "hermidalc@pitt.edu"
__license__ = "MIT"

import re

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq") or snakemake.input.get("fq1")
assert fq1 is not None, "input: fq/fq1 is a required input parameter"
in_fqs = f"--in1 {fq1}"

fq2 = snakemake.input.get("fq2")
if fq2:
    in_fqs += f" --in2 {fq2}"

trim_fq1 = snakemake.output.get("trim_fq1")
if trim_fq1:
    out_fqs = f"--out1 {trim_fq1}"
    trim_fq2 = snakemake.output.get("trim_fq2")
    if trim_fq2:
        out_fqs += f" --out2 {trim_fq2}"
        trim_up1 = snakemake.output.get("trim_up1")
        if trim_up1:
            out_fqs += f" --unpaired1 {trim_up1}"
        trim_up2 = snakemake.output.get("trim_up2")
        if trim_up2:
            out_fqs += f" --unpaired2 {trim_up2}"
        merged = snakemake.output.get("merged")
        if merged:
            if not re.search(r"--merge\b", extra):
                raise ValueError(
                    "output.merged specified but '--merge' option missing from params.extra"
                )
            out_fqs += f" --merged_out {merged}"
else:
    out_fqs = ""

failed = snakemake.output.get("failed")
if failed:
    out_fqs += f" --failed_out {failed}"

adapters = snakemake.params.get("adapters", "")

extra = snakemake.params.get("extra", "")

html = snakemake.output.get("html")
assert html is not None, "output: html is a required output parameter"
json = snakemake.output.get("json")
assert json is not None, "output: json is a required output parameter"

run_id = snakemake.params.get("run_id")
title = f"{run_id} fastp report" if run_id else "fastp report"

shell(
    "fastp"
    " --thread {snakemake.threads}"
    " {in_fqs}"
    " {out_fqs}"
    " {adapters}"
    " {extra}"
    " --html {html}"
    " --json {json}"
    " --report_title {title:q}"
    " {log}"
)
