# ===============================================
# Module 6: FEELnc – Coding potential & classification
# ===============================================

# ---------------------------
# FEELnc filter: remove transcripts overlapping protein‑coding regions
# ---------------------------
rule feelnc_filter:
    input:
        candidates = config["candidate_gtf"],
        ref_gtf    = config["genome_gtf"],
        genome_fa  = config["genome_fasta"]
    output:
        filtered = config["feelnc_filter_out"]
    benchmark:
        "benchmarks/feelnc_filter.txt"
    conda:
        "../envs/feelnc_full.yaml"
    shell:
        r"""
        export FEELNCPATH={workflow.basedir}/tools/FEELnc
        export PERL5LIB=$FEELNCPATH/lib
        export PATH=$FEELNCPATH/scripts:$FEELNCPATH/utils:$FEELNCPATH/bin/LINUX:$PATH

        mkdir -p results/feelnc
        perl $FEELNCPATH/scripts/FEELnc_filter.pl \
          -i {input.candidates} \
          -a {input.ref_gtf} \
          -b transcript_biotype=protein_coding \
          > {output.filtered}
        """

# ---------------------------
# FEELnc coding potential (random forest model)
# ---------------------------
rule feelnc_codpot:
    input:
        filtered    = config["feelnc_filter_out"],
        genome      = config["genome_fasta"],
        train_mrna  = config["train_mrna_fasta"],
        train_lncrna= config["train_lncrna_fasta"]
    output:
        lncrna_gtf = config["feelnc_lncrna_out"]
    benchmark:
        "benchmarks/feelnc_codpot.txt"
    conda:
        "../envs/feelnc_full.yaml"
    shell:
        r"""
        export FEELNCPATH={workflow.basedir}/tools/FEELnc
        export PERL5LIB=$FEELNCPATH/lib
        export PATH=$FEELNCPATH/scripts:$FEELNCPATH/utils:$FEELNCPATH/bin/LINUX:$PATH

        mkdir -p results/feelnc
        perl $FEELNCPATH/scripts/FEELnc_codpot.pl \
          -i {input.filtered} \
          -a {input.train_mrna} \
          -l {input.train_lncrna} \
          -g {input.genome} \
          --outdir results/feelnc
        """

# ---------------------------
# FEELnc classifier: genomic context classification
# ---------------------------
rule feelnc_classifier:
    input:
        gtf = config["feelnc_lncrna_out"],
        annotation = config["genome_gtf"]
    output:
        classes = config["feelnc_classifier_out"]
    benchmark:
        "benchmarks/feelnc_classifier.txt"
    conda:
        "../envs/feelnc_full.yaml"
    shell:
        r"""
        export FEELNCPATH={workflow.basedir}/tools/FEELnc
        export PERL5LIB=$FEELNCPATH/lib
        export PATH=$FEELNCPATH/scripts:$FEELNCPATH/utils:$FEELNCPATH/bin/LINUX:$PATH
        perl $FEELNCPATH/scripts/FEELnc_classifier.pl \
          -i {input.gtf} \
          -a {input.annotation} \
          > {output.classes}
        """
