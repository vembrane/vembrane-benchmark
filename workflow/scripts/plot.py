import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots

input_tsv = snakemake.input.data
scenarios = snakemake.params.scenarios

df = pd.read_csv(input_tsv, sep="\t")
df["seconds_per_record"] = 1.0 / df["records_per_second"]
df = df[df["scenario"].isin(scenarios)]

non_plain_vcf = df.query(
    "(tool == 'vembrane' | tool == 'bcftools' | tool == 'slivar') & filetype != 'vcf_v' & filetype != 'vcf_z'"
)
df = df.drop(non_plain_vcf.index)

bcf_tools = ["vembrane", "bcftools", "slivar"]
for tool in bcf_tools:
    non_plain_vcf.loc[
        non_plain_vcf.query(f"tool == '{tool}'").index, "tool"
    ] = f"{tool} (bcf)"
df = df.append(non_plain_vcf)
df["n_samples"] = df["sample"].str.count("_") + 1
tools = set(df["tool"])

# paul tol's bright qualitative
palette = [
    "#4477AA",  # bcftools (bcf);  blue
    "#7198bd",  # bcftools      ;  blue, s-20, v+7
    "#EE6677",  # vembrane (bcf);  red
    "#ffa1ac",  # vembrane      ;  red,  s-20, v+7
    "#228833",  # SnpSift       ;  green
    "#CCBB44",  # bio-vcf       ;  yellow
    "#66CCEE",  # slivar (bcf)  ;  cyan
    "#a1e7ff",  # slivar        ;  cyan,  s-20, v+7
    "#AA3377",  # filter_vep    ;  purple
]

fig = px.violin(
    df,
    y="records_per_second",
    color="tool",
    # box=True,
    # points=False,
    facet_col="scenario",
    # facet_row="sample",
    template="plotly_white",
    color_discrete_sequence=palette,
    category_orders={
        "tool": [
            "bcftools (bcf)",
            "bcftools",
            "vembrane (bcf)",
            "vembrane",
            "SnpSift",
            "bio-vcf",
            "slivar (bcf)",
            "slivar",
            "filter_vep",
        ],
        "scenario": scenarios,
    },
    labels={
        "seconds_per_record": "seconds per record",
        "tool": "",
        "records_per_second": "records per second",
    },
    log_y=True,
)

for tool in tools:
    if tool in bcf_tools:
        fig.update_traces(
            ## sadly plotly.graph_object.Violin traces have no "dash" property on the line markers.
            # patch=dict(line=dict(dash="dot")),
            legendgroup=tool,
            selector=dict(name=f"{tool} (bcf)"),
        )
    else:
        fig.update_traces(
            legendgroup=tool,
            selector=dict(name=tool),
        )

fig.update_traces(spanmode="hard", selector=dict(type="violin"))
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
fig.for_each_trace(
    lambda t: t.update(
        scalemode="width", scalegroup=f"{id(t)}", meanline=dict(visible=False)
    )
)
fig.update_xaxes(zeroline=True, showline=True, showgrid=False)
fig.update_layout(
    violingap=0.0,
    violingroupgap=0.05,
    boxgap=0.0,
    boxgroupgap=0.1,
    margin=dict(l=0, r=0, t=20, b=0),
)


# from https://github.com/plotly/plotly.py/issues/1930#issuecomment-645079975
def modify_facet_spacing(fig, ncols, hs, vs):
    import re, math, numpy as np

    pattern = r"[xy]axis(\d+)?"

    fct_titles = fig.layout["annotations"]
    nrows = int(
        math.ceil(len(fct_titles) / ncols)
    )  # this needs revision for cases where plotly express has removed redundant rows (if for eg the bottom row or two are empty of charts)

    sp = make_subplots(
        rows=nrows, cols=ncols, horizontal_spacing=hs, vertical_spacing=vs
    )

    # the subplots arange plots differently than faceted plots
    # subplots: top left, facets: bottom left
    # pdb.set_trace()
    spt_rng = np.arange(1, (nrows * ncols) + 1)
    fct_rng = np.flip(spt_rng.reshape(nrows, ncols), axis=0).flatten()
    repl = dict(zip(fct_rng, spt_rng))

    empty_cells = (nrows * ncols) - len(fct_titles)

    for el in fig.layout:

        if re.match(pattern, el):
            axis_fct = el
            axis_fct_n = int(re.match(pattern, axis_fct).groups(0)[0])

            if "yaxis" in axis_fct:
                if axis_fct_n == 0:
                    # pdb.set_trace()
                    axis_fct_n = 1
                    axis_fct += str(axis_fct_n)
                axis_spt = axis_fct.replace(str(axis_fct_n), str(repl[axis_fct_n]))
            else:
                axis_spt = axis_fct

            spt_index = int(re.match(pattern, axis_spt).groups(0)[0]) - 1
            fct_index = axis_fct_n - 1

            if axis_fct.endswith("axis1"):
                axis_fct = axis_fct[:-1]
            if axis_spt.endswith("axis1"):
                axis_spt = axis_spt[:-1]

            fig.layout[axis_fct].update(domain=sp.layout[axis_spt].domain)

            if "yaxis" in axis_fct:
                # facets are nrows*ncols in size behind the scenes, but
                # fct_titles is only ncharts long (eg if 1 plot on last row) -
                # the rest of the plots on that row are empty
                if fct_index - empty_cells >= 0:
                    fct_titles[fct_index - empty_cells]["y"] = sp.layout[
                        axis_spt
                    ].domain[1]

    fig.layout["annotations"] = fct_titles


modify_facet_spacing(fig, len(scenarios), 0.015, 0)

for img_path in snakemake.output.images:
    fig.write_image(img_path, width=1000, height=400)
fig.write_html(snakemake.output.html)
