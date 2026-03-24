rule stringtie_pass1:
    input:
        bam = "results/aligned_bam/{sample}.sorted.bam",
        bai = "results/aligned_bam/{sample}.sorted.bam.bai",
        gtf = config["genome_gtf"]
    output:
        "results/assembly/{sample}.gtf"
    benchmark:
        "benchmarks/stringtie_pass1/{sample}.txt"
    threads: 4
    conda: "../envs/assembly.yaml"
    shell:
        "stringtie {input.bam} -G {input.gtf} -o {output} -p {threads}"
