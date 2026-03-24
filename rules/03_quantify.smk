# ===============================================
# Module 3: Quantification (Count Matrices)
# ===============================================

import pandas as pd

samples = pd.read_table(config["samples"], dtype=str)

rule prepDE:
    """Generate gene and transcript count matrices using prepDE.py in directory mode."""
    input:
        # The rule depends on all sample GTFs from stringtie_pass3
        expand("results/quant/{sample}/{sample}.gtf", sample=samples["sample"])
    output:
        gene_counts = config["gene_counts"],
        transcript_counts = config["transcript_counts"]
    params:
        script = "scripts/prepDE.py",
        quant_dir = "results/quant"
    conda:
        "../envs/prepDE.yaml"
    run:
        # Run prepDE.py pointing to the quant directory
        shell("python {params.script} -i {params.quant_dir}")
        # Move the generated files to the expected output locations
        shell("mv gene_count_matrix.csv {output.gene_counts}")
        shell("mv transcript_count_matrix.csv {output.transcript_counts}")
