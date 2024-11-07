from os.path import basename, join


def get_fq(wildcards, data_type="raw"):
    u = UNIT_DF.loc[
        (
            (wildcards.sample, wildcards.unit)
            if hasattr(wildcards, "unit")
            else wildcards.sample
        ),
        ["fq1", "fq2"],
    ].map(
        lambda x: (
            join(
                (
                    HOST_FILTER_RESULTS_DIR
                    if data_type == "filtered"
                    else (
                        TRIMMED_RESULTS_DIR
                        if config["trim"]["activate"] and data_type == "trimmed"
                        else join(FASTQ_DATA_DIR, x)
                    )
                ),
                basename(x),
            )
        )
    )
    return (
        {"reads": [u.fq1, u.fq2]}
        if data_type == "trimmed"
        else {"fq1": u.fq1, "fq2": u.fq2}
    )
