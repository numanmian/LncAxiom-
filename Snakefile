import pandas as pd

# Load configuration
configfile: "config/config.yaml"

# Include modules based on config and input method
if config.get("run_download", False) and config.get("input_method", "sra") == "sra":
    include: "rules/00_download_sra.smk"
elif config.get("run_download", False) and config.get("input_method") == "local":
    print("Warning: run_download is true but input_method is 'local'. Download rule will be skipped.")

if config.get("run_preprocessing", False):
    include: "rules/01_preprocess.smk"

if config.get("run_assembly", False):
    include: "rules/02_assemble.smk"

if config.get("run_quantification", False):
    include: "rules/03_quantify.smk"

if config.get("run_lncrna_prediction", False):
    include: "rules/04_lncrna_prediction.smk"

# Optional FEELnc module
if (
    config.get("run_feelnc_filter", False)
    or config.get("run_feelnc_codpot", False)
    or config.get("run_feelnc_classifier", False)
    or config.get("run_feelnc_intersection", False)
):
    include: "rules/06_feelnc.smk"

# psRobot eTM prediction module
if config.get("run_psrobot_etm", False):
    include: "rules/08_psrobot_etm.smk"

# BLAST novelty assessment module
if config.get("run_blast_novelty", False):
    include: "rules/09_blast_novelty.smk"

# Functional enrichment module (Pearson correlation + GO/KEGG)
if config.get("run_enrichment", False):
    include: "rules/11_functional_enrichment.smk"

# DESeq2 differential expression module
include: "rules/deseq2.smk"

# ----------------------------------------------------
# Default target: run all enabled modules
# ----------------------------------------------------
rule all:
    input:
        # Preprocessing outputs
        expand("results/multiqc/multiqc_report.html") if config.get("run_preprocessing", False) else [],
        expand("results/trimmed_fastq/{sample}.fastq.gz",
            sample=pd.read_table(config["samples"], dtype=str)["sample"],
        ) if config.get("run_preprocessing", False) else [],

        # Assembly outputs
        expand("results/assembly/{sample}.gtf",
            sample=pd.read_table(config["samples"], dtype=str)["sample"],
        ) if config.get("run_assembly", False) else [],
        expand("results/gffcompare/merged.{ext}",
            ext=["annotated.gtf", "stats", "tracking", "loci"],
        ) if config.get("run_assembly", False) else [],

        # Quantification outputs
        config["gene_counts"]
        if config.get("run_quantification", False)
        else [],
        config["transcript_counts"]
        if config.get("run_quantification", False)
        else [],

        # lncRNA prediction outputs (gffcompare-based)
        config.get("candidate_gtf", "results/lncrna_candidates.gtf")
        if config.get("run_lncrna_prediction", False)
        else [],

        # FEELnc outputs (only if enabled in config)
        config["feelnc_filter_out"]
        if config.get("run_feelnc_filter", False)
        else [],
        config["feelnc_lncrna_out"]
        if config.get("run_feelnc_codpot", False)
        else [],
        config["feelnc_classifier_out"]
        if config.get("run_feelnc_classifier", False)
        else [],
        config["feelnc_intersection_out"]
        if config.get("run_feelnc_intersection", False)
        else [],

        # DESeq2 outputs (always generated)
        expand("results/{level}_deseq2_results.csv", level=["gene", "transcript"]),
        expand("plots/{level}_volcano.png", level=["gene", "transcript"]),

        # psRobot eTM outputs
        config["psrobot_output_dir"] + "/etm_predictions.tsv" if config.get("run_psrobot_etm", False) else [],

        # BLAST novelty assessment outputs
        config["blast_output_dir"] + "/novelty_summary.tsv" if config.get("run_blast_novelty", False) else [],

        # Functional enrichment outputs (Pearson correlation + GO/KEGG)
        expand("results/enrichment/{file}",
               file=["significant_correlations.csv",
                     "enrichR_GO_Biological_Process_2023.csv",
                     "enrichR_Reactome_Pathways_2024.csv"]) if config.get("run_enrichment", False) else []
