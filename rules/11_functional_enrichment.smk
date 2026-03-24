# ===============================================
# Module 11: Functional Enrichment (Pearson correlation + GO/KEGG)
# ===============================================

rule functional_enrichment:
    """Pearson correlation between DE lncRNAs and mRNAs, then enrichment of co-expressed mRNAs."""
    input:
        gene_counts   = config["gene_counts"],
        deseq2_results= "results/gene_deseq2_results.csv",
        lncrna_gtf    = config["feelnc_lncrna_out"]   # path to filter.gtf.lncRNA.gtf
    output:
        corr_table    = "results/enrichment/significant_correlations.csv",
        go_enrich     = "results/enrichment/enrichR_GO_Biological_Process_2023.csv",
        reactome_enrich= "results/enrichment/enrichR_Reactome_Pathways_2024.csv"
    params:
        cor_threshold = config.get("enrichment_cor_threshold", 0.8),
        p_threshold   = config.get("enrichment_p_threshold", 0.01),
        databases     = config.get("enrichment_databases",
                                   ["GO_Biological_Process_2023",
                                    "Reactome_Pathways_2024"])
    conda:
        "../envs/enrichment.yaml"
    script:
        "../scripts/functional_enrichment.R"
