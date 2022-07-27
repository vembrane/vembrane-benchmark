rule plot:
    input:
        data="results/benchmark.tsv",
    output:
        html="results/plot/benchmark.html",
        images=expand("results/plot/benchmark.{ext}", ext=["pdf", "svg"]),
    params:
        scenarios=SCENARIOS,
    log:
        "logs/plot/plot.log",
    conda:
        "../envs/plotly.yml"
    script:
        "../scripts/plot.py"


rule aggregate_benchmarks:
    input:
        files=expand(
            "benchmarks/aux/{item.annotation}/{item.tool}/{item.scenario}/{item.filetype[0]}/{item.sample}_annotated_{item.annotation}.{item.filetype[1]}.benchmark.tsv",
            item=VALID_WILDCARD_COMBINATIONS.query("sample != 'HG002'").itertuples(),
        ),
    output:
        tsv="results/benchmark.tsv",
    log:
        "logs/plot/aggregate_benchmarks.log",
    run:
        import pandas as pd

        df = pd.concat(pd.read_csv(f, sep="\t", index_col=None) for f in input.files)
        df.to_csv(output.tsv, sep="\t", index=False)


rule add_benchmark_info:
    input:
        f="benchmarks/filtered/{annotation}/{tool}/{scenario}/{filetype}/{sample}_annotated_{annotation}.{extension}.benchmark.txt",
        counts_total="results/counts/{sample}.txt",
        counts_out="results/counts_out/{annotation}/{tool}/{scenario}/{filetype}/{sample}_annotated_{annotation}.{extension}.txt",
    output:
        f="benchmarks/aux/{annotation}/{tool}/{scenario}/{filetype}/{sample}_annotated_{annotation}.{extension}.benchmark.tsv",
    conda:
        "../envs/plain.yaml"
    log:
        "logs/plot/add_benchmark_info/{annotation}.{tool}.{scenario}.{filetype}.{sample}.{annotation}.{extension}.log",
    threads: 1
    script:
        "../scripts/add_benchmark_info.py"
