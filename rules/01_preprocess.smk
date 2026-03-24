# ===============================================
# Module 1: Preprocessing (Single-End Only)
# ===============================================

import pandas as pd

# Read sample table
samples = pd.read_table(config["samples"], dtype=str)

# Helper to get raw FASTQ path based on input method
def get_raw_fastq(wildcards):
    """Return path to raw FASTQ:
       - SRA mode: data/raw_fastq/{sra_id}.fastq.gz
       - Local mode: {local_fastq_dir}/{sample}.fastq.gz
    """
    sample = wildcards.sample
    if config.get("input_method", "sra") == "sra":
        # Use SRA ID from samples table
        sra_id = samples[samples["sample"] == sample]["sra"].values[0]
        return f"data/raw_fastq/{sra_id}.fastq.gz"
    else:
        # Local mode: use local_fastq_dir and sample name
        local_dir = config.get("local_fastq_dir", "data/raw_fastq")
        return f"{local_dir}/{sample}.fastq.gz"

# ---------------------------
# Rule: FastQC (raw)
# ---------------------------
rule fastqc_raw:
    input:
        fastq = get_raw_fastq
    output:
        html = "results/qc/fastqc_raw/{sample}_fastqc.html",
        zip = "results/qc/fastqc_raw/{sample}_fastqc.zip"
    benchmark:
        "benchmarks/fastqc_raw/{sample}.txt"
    params:
        threads = config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "fastqc -t {params.threads} -o results/qc/fastqc_raw {input.fastq}"

# ---------------------------
# Rule: fastp trimming
# ---------------------------
rule fastp_trim:
    input:
        r1 = get_raw_fastq
    output:
        r1_out = "results/trimmed_fastq/{sample}.fastq.gz",
        json = "results/qc/fastp/{sample}.json",
        html = "results/qc/fastp/{sample}.html"
    benchmark:
        "benchmarks/fastp_trim/{sample}.txt"
    params:
        args = config.get("fastp_args", ""),
        threads = config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "fastp {params.args} -i {input.r1} -o {output.r1_out} "
        "-j {output.json} -h {output.html} -w {params.threads}"

# ---------------------------
# Rule: FastQC (trimmed)
# ---------------------------
rule fastqc_trimmed:
    input:
        r1 = "results/trimmed_fastq/{sample}.fastq.gz"
    output:
        html = "results/qc/fastqc_trimmed/{sample}_fastqc.html",
        zip = "results/qc/fastqc_trimmed/{sample}_fastqc.zip"
    benchmark:
        "benchmarks/fastqc_trimmed/{sample}.txt"
    params:
        threads = config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "fastqc -t {params.threads} -o results/qc/fastqc_trimmed {input.r1}"

# ---------------------------
# Rule: HISAT2 index (build if not exists)
# ---------------------------
rule hisat2_index:
    input:
        fasta = config["genome_fasta"]
    output:
        touch("results/indexes/TAIR10.done")
    benchmark:
        "benchmarks/hisat2_index.txt"
    params:
        index_base = config["hisat2_index"],
        threads = config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "hisat2-build -p {params.threads} {input.fasta} {params.index_base} && touch {output}"

# ---------------------------
# Rule: HISAT2 alignment (single-end)
# ---------------------------
rule hisat2_align:
    input:
        r1 = "results/trimmed_fastq/{sample}.fastq.gz",
        idx_done = "results/indexes/TAIR10.done"
    output:
        sam = "results/aligned_bam/{sample}.sam"
    benchmark:
        "benchmarks/hisat2_align/{sample}.txt"
    params:
        args = config.get("hisat2_args", ""),
        threads = config["threads"],
        index_base = config["hisat2_index"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "hisat2 -p {params.threads} -x {params.index_base} -U {input.r1} {params.args} -S {output.sam}"

# ---------------------------
# Rule: samtools sort
# ---------------------------
rule samtools_sort:
    input:
        "results/aligned_bam/{sample}.sam"
    output:
        bam = "results/aligned_bam/{sample}.sorted.bam"
    benchmark:
        "benchmarks/samtools_sort/{sample}.txt"
    params:
        threads = config["threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "samtools sort -@ {params.threads} -o {output.bam} {input}"

# ---------------------------
# Rule: samtools index
# ---------------------------
rule samtools_index:
    input:
        "results/aligned_bam/{sample}.sorted.bam"
    output:
        "results/aligned_bam/{sample}.sorted.bam.bai"
    benchmark:
        "benchmarks/samtools_index/{sample}.txt"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "samtools index {input}"

# ---------------------------
# Rule: MultiQC (all QC)
# ---------------------------
rule multiqc:
    input:
        expand("results/qc/fastqc_raw/{sample}_fastqc.zip", sample=samples["sample"]),
        expand("results/qc/fastqc_trimmed/{sample}_fastqc.zip", sample=samples["sample"]),
        expand("results/qc/fastp/{sample}.json", sample=samples["sample"])
    output:
        html = "results/multiqc/multiqc_report.html"
    benchmark:
        "benchmarks/multiqc.txt"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        "multiqc --force -o results/multiqc results/qc"
