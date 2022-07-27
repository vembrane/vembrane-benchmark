import pandas as pd


counts_total = 0
with open(snakemake.input.counts_total, "rt") as countfile:
    counts_total = int(countfile.read())

counts_out = 0
with open(snakemake.input.counts_out, "rt") as countfile:
    counts_out = int(countfile.read())

df = pd.read_csv(
    snakemake.input.f,
    sep="\t",
    index_col=None,
)

wc = snakemake.wildcards
df["scenario"] = wc.scenario
df["sample"] = wc.sample
df["filetype"] = wc.filetype
df["tool"] = wc.tool
df["annotation"] = wc.annotation
df["count"] = counts_total
df["records_written"] = counts_out
df["records_per_second"] = counts_total / df["s"]
df.to_csv(snakemake.output.f, sep="\t", index=False)