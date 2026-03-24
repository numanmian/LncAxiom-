# ===============================================
# Module: psRobot eTM prediction
# ===============================================

# Rule to convert final lncRNA GTF to FASTA
rule gtf_to_fasta_lncrna:
    input:
        gtf = config["feelnc_final_gtf"],
        genome = config["genome_fasta"]
    output:
        fasta = config["feelnc_final_gtf"].replace(".gtf", ".fa")
    conda:
        "../envs/gffread.yaml"   # we need an environment with gffread
    shell:
        "gffread -w {output.fasta} -g {input.genome} {input.gtf}"

# psRobot eTM prediction
rule psrobot_etm:
    input:
        lncrna_fasta = "results/feelnc/candidates.fa",
        mirna_fasta = config["psrobot_mirna_fasta"]
    output:
        etm_tsv = config["psrobot_output_dir"] + "/etm_predictions.tsv"
    params:
        score_thresh = config.get("psrobot_score_thresh", 2.5),
        seed_start = 2,
        seed_end = 8,
        forbidden_bulge_start = 9,
        forbidden_bulge_end = 12,
        max_bulge_size = config.get("psrobot_max_bulge_size", 3),
        max_mismatches_gu = config.get("psrobot_max_mismatches_gu", 3),
        threads = config.get("threads", 8),
        psrobot_bin = config.get("psrobot_bin", "psRobot_tar"),
        script = "scripts/run_psrobot_etm.py"
    conda:
        "../envs/psrobot.yaml"   # we will create this environment next
    shell:
        "python {params.script} "
        "--lncrna_fasta {input.lncrna_fasta} "
        "--mirna_fasta {input.mirna_fasta} "
        "--output {output.etm_tsv} "
        "--score_thresh {params.score_thresh} "
        "--seed_start {params.seed_start} "
        "--seed_end {params.seed_end} "
        "--forbidden_bulge_start {params.forbidden_bulge_start} "
        "--forbidden_bulge_end {params.forbidden_bulge_end} "
        "--max_bulge_size {params.max_bulge_size} "
        "--max_mismatches_gu {params.max_mismatches_gu} "
        "--threads {params.threads} "
        "--psrobot_bin {params.psrobot_bin}"
