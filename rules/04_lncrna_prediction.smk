# ===============================================
# Module 4: lncRNA Prediction
# ===============================================

import pandas as pd

samples = pd.read_table(config["samples"], dtype=str)

# Rule: filter transcripts by class code (i,u,x) from annotated GTF
rule filter_classcode:
    """Extract transcripts with class codes i, u, x from gffcompare annotated GTF."""
    input:
        annotated_gtf = "results/gffcompare/merged.annotated.gtf"
    output:
        candidate_gtf = config.get("candidate_gtf", "results/lncrna_candidates.gtf")
    params:
        script = "scripts/filter_by_classcode.py"
    run:
        shell("python {params.script} {input.annotated_gtf} {output}")
