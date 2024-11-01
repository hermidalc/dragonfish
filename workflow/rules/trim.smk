rule fastp_trimmed_fastq:
    input:
        unpack(lambda wc: get_fq(wc, trimmed=False)),
    params:
        run_id=RUN_ID_WILDCARD_STR,
        extra=config["trim"]["fastp"]["extra"],
    output:
        trim_fq1=TRIMMED_FASTQ1_FILE,
        trim_fq2=TRIMMED_FASTQ2_FILE,
        trim_up1=TRIMMED_UNPAIR1_FILE,
        trim_up2=TRIMMED_UNPAIR2_FILE,
        failed=FAILED_READS_FILE,
        html=FASTP_HTML_REPORT_FILE,
        json=FASTP_JSON_REPORT_FILE,
    threads: config["trim"]["fastp"]["threads"]
    log:
        FASTP_LOG,
    wrapper:
        FASTP_WRAPPER
