# ===============================================
# Module 0: SRA Download and FASTQ Conversion
# ===============================================

import os

# Rule: prefetch_sra
rule prefetch_sra:
    """Download SRA file using prefetch."""
    output:
        sra = "data/sra_cache/{sra}/{sra}.sra"
    benchmark:
        "benchmarks/prefetch_sra/{sra}.txt"
    params:
        outdir = "data/sra_cache"
    log:
        "logs/prefetch/{sra}.log"
    conda:
        "../envs/sra_tools.yaml"
    run:
        shell("prefetch {wildcards.sra} -O {params.outdir} > {log} 2>&1")

# Rule: fasterq_dump (single-end)
rule fasterq_dump:
    """Convert .sra to FASTQ and compress."""
    input:
        sra = "data/sra_cache/{sra}/{sra}.sra"
    output:
        fastq = "data/raw_fastq/{sra}.fastq.gz"
    benchmark:
        "benchmarks/fasterq_dump/{sra}.txt"
    params:
        outdir = "data/raw_fastq",
        threads = config["threads"]
    log:
        "logs/fasterq_dump/{sra}.log"
    conda:
        "../envs/sra_tools.yaml"
    run:
        shell("mkdir -p {params.outdir}")
        shell("fasterq-dump {input.sra} -O {params.outdir} -e {params.threads} > {log} 2>&1")
        shell("pigz -p {params.threads} {params.outdir}/{wildcards.sra}.fastq")

# Rule: cleanup_sra (optional)
rule cleanup_sra:
    """Remove .sra file after successful conversion to save space."""
    input:
        fastq = "data/raw_fastq/{sra}.fastq.gz"
    output:
        touch("data/sra_cache/{sra}/.deleted")
    benchmark:
        "benchmarks/cleanup_sra/{sra}.txt"
    run:
        sra_file = f"data/sra_cache/{wildcards.sra}/{wildcards.sra}.sra"
        if os.path.exists(sra_file):
            os.remove(sra_file)
        try:
            os.rmdir(os.path.dirname(sra_file))
        except OSError:
            pass
