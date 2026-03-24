# ===============================================
# Module 2: Transcript Assembly (StringTie 3-pass)
# ===============================================

import pandas as pd

samples = pd.read_table(config["samples"], dtype=str)

# ---------------------------
# PASS 1: Reference-guided assembly per sample
# ---------------------------
rule stringtie_pass1:
    """Assemble transcripts per sample (allowing novel transcripts)."""
    input:
        bam = "results/aligned_bam/{sample}.sorted.bam",
        bai = "results/aligned_bam/{sample}.sorted.bam.bai",
        gtf = config["genome_gtf"]
    output:
        gtf = "results/assembly/{sample}.gtf"
    params:
        threads = config["threads"],
        extra   = config.get("stringtie_args", "")
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        stringtie {input.bam} \
          -p {params.threads} \
          -G {input.gtf} \
          -o {output.gtf} \
          {params.extra}
        """

# ---------------------------
# PASS 2: Merge all assemblies into a unified transcriptome
# ---------------------------
rule stringtie_merge:
    """Merge all sample assemblies with reference annotation."""
    input:
        assemblies = expand("results/assembly/{sample}.gtf", sample=samples["sample"]),
        ref_gtf    = config["genome_gtf"]
    output:
        merged = config["merged_gtf"]
    params:
        threads = config["threads"],
        extra   = config.get("stringtie_args", "")
    conda:
        "../envs/stringtie.yaml"
    run:
        list_file = "results/assembly/assembly_list.txt"
        with open(list_file, "w") as f:
            for gtf in input.assemblies:
                f.write(gtf + "\n")
        shell(
            "stringtie --merge -p {threads} -G {ref_gtf} "
            "-o {merged} {extra} {list_file}"
            .format(
                threads=params.threads,
                ref_gtf=input.ref_gtf,
                merged=output.merged,
                extra=params.extra,
                list_file=list_file,
            )
        )

# ---------------------------
# PASS 3: Quantify against merged transcriptome
# ---------------------------
rule stringtie_pass3:
    """Quantify using merged GTF (no novel assembly)."""
    input:
        bam        = "results/aligned_bam/{sample}.sorted.bam",
        bai        = "results/aligned_bam/{sample}.sorted.bam.bai",
        merged_gtf = config["merged_gtf"]
    output:
        quant_dir = directory("results/quant/{sample}"),
        gtf       = "results/quant/{sample}/{sample}.gtf"
    params:
        threads = config["threads"],
        extra   = "-e -B"  # estimate only + Ballgown tables
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        mkdir -p results/quant/{wildcards.sample}
        stringtie {input.bam} \
          -p {params.threads} \
          -G {input.merged_gtf} \
          -o {output.gtf} \
          {params.extra}
        """

# ---------------------------
# Rule: gffcompare (compare merged transcripts to reference)
# ---------------------------
rule gffcompare:
    """Compare merged transcripts to reference annotation and assign class codes."""
    input:
        merged_gtf = config["merged_gtf"],
        ref_gtf    = config["genome_gtf"]
    output:
        annotated = "results/gffcompare/merged.annotated.gtf",
        stats     = "results/gffcompare/merged.stats",
        tracking  = "results/gffcompare/merged.tracking",
        loci      = "results/gffcompare/merged.loci"
    params:
        outprefix = "results/gffcompare/merged"
    conda:
        "../envs/stringtie.yaml"
    shell:
        """
        mkdir -p results/gffcompare
        gffcompare \
          -r {input.ref_gtf} \
          -o {params.outprefix} \
          {input.merged_gtf}
        """
