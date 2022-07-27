from snakemake.io import repeat


rule filter_directly:
    input:
        "results/{annotation}/{filetype}/{file}",
    output:
        "results/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.vcf",
    benchmark:
        repeat(
            "benchmarks/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.benchmark.txt",
            BENCHMARK_REPEATS,
        )
    log:
        "logs/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.log.txt",
    params:
        cmd=lambda wc: CMD_TEMPLATES["direct"][wc.tool].substitute(
            FILTER=escape_expression(
                CONFIGS[wc.tool]["annotations"][wc.annotation][wc.scenario]
            )
        ),
    conda:
        "../envs/{tool}.yml"
    threads: 1
    shell:
        """
        {params.cmd} {input} > {output} 2> {log}
        """
