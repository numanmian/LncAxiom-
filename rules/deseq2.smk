rule deseq2_with_volcano:
    input:
        counts   = lambda wc: config["counts"][wc.level],
        metadata = config["metadata"]
    output:
        csv     = "results/{level}_deseq2_results.csv",
        volcano = "plots/{level}_volcano.png"
    params:
        design      = config["design"],
        factor      = config["contrast"]["factor"],
        case        = config["contrast"]["case"],
        control     = config["contrast"]["control"],
        padj_cutoff = config["volcano"]["padj_cutoff"],
        lfc_cutoff  = config["volcano"]["lfc_cutoff"],
        level       = "{level}"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/run_deseq2_and_volcano.R"
