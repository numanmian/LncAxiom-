# ===============================================
# Module 3: Merge assemblies & quantification
# ===============================================

import pandas as pd

# Read sample table to get list of samples
samples = pd.read_table(config["samples"], dtype=str)

# ---------------------------
# Rule: Merge sample assemblies
# ---------------------------
rule stringtie_merge:
    input:
        gtf_list = expand("results/assembly/{sample}.gtf", sample=samples["sample"]),
        ref_gtf = config["genome_gtf"]
    output:
        merged = config["merged_gtf"]
    benchmark:
        "benchmarks/stringtie_merge.txt"
    conda:
        "../envs/assembly.yaml"
    shell:
        "stringtie --merge -G {input.ref_gtf} -o {output.merged} {input.gtf_list}"

# ---------------------------
# Rule: StringTie pass3 (quantification)
# ---------------------------
rule stringtie_pass3:
    input:
        bam = "results/aligned_bam/{sample}.sorted.bam",
        bai = "results/aligned_bam/{sample}.sorted.bam.bai",
        merge_gtf = config["merged_gtf"]
    output:
        dir = directory("results/quant/{sample}"),
        gtf = "results/quant/{sample}/{sample}.gtf"
    benchmark:
        "benchmarks/stringtie_pass3/{sample}.txt"
    threads: 4
    conda:
        "../envs/assembly.yaml"
    shell:
        "stringtie {input.bam} -G {input.merge_gtf} -o {output.gtf} -p {threads} -e -B"

# ---------------------------
# Rule: Generate count matrices (prepDE)
# ---------------------------
rule prepDE:
    input:
        expand("results/quant/{sample}/{sample}.gtf", sample=samples["sample"])
    output:
        gene = config["gene_counts"],
        transcript = config["transcript_counts"]
    benchmark:
        "benchmarks/prepDE.txt"
    conda:
        "../envs/assembly.yaml"
    shell:
        "prepDE.py -i results/quant -g {output.gene} -t {output.transcript}"
