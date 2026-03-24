# ===============================================
# Module: BLAST novelty assessment
# ===============================================

# Rule to run BLAST against each database
rule blast_against_db:
    input:
        query = "results/feelnc/candidates.fa",
        db_index = lambda wildcards: config["blast_databases"][wildcards.db] + ".ndb"
    output:
        "results/blast_novelty/{db}.blast"
    params:
        db_base = lambda wildcards: config["blast_databases"][wildcards.db],
        evalue = config["blast_evalue"],
        outfmt = config["blast_outfmt"],
        threads = config["blast_threads"]
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn -query {input.query} -db {params.db_base} -out {output} "
        "-outfmt '{params.outfmt}' -evalue {params.evalue} -num_threads {params.threads}"

# Rule to combine all BLAST results into a summary TSV
rule summarize_blast:
    input:
        blast_files = expand("results/blast_novelty/{db}.blast", db=config["blast_databases"].keys())
    output:
        summary = "results/blast_novelty/novelty_summary.tsv"
    params:
        dbs = list(config["blast_databases"].keys()),
        script = "scripts/summarize_blast.py"
    conda:
        "../envs/blast.yaml"
    shell:
        "python {params.script} --blast_dir results/blast_novelty --dbs {params.dbs} --output {output.summary}"
