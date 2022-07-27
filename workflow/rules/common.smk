import pandas as pd

from string import Template
from itertools import product, combinations

CONFIGS = config["benchmark"]["tool"]
TOOLS = list(CONFIGS.keys())

CMD_TEMPLATES = {
    "direct": {tool: Template(data["invocation"]) for tool, data in CONFIGS.items()}
}

SAMPLE_NUMS = [1, 2, 3, 4]

SAMPLES = [
    "_".join(c)
    for k in range(1, len(SAMPLE_NUMS) + 1)
    for c in combinations([f"HG00{i}" for i in SAMPLE_NUMS], k)
]

FILETYPES = [("bcf_b", "bcf"), ("bcf_u", "bcf"), ("vcf_v", "vcf"), ("vcf_z", "vcf.gz")]

VALID_WILDCARD_COMBINATIONS = []
SCENARIOS = []

for tool, data in CONFIGS.items():
    annotations = data["annotations"]
    for annotation, scenarios in annotations.items():
        for scenario, expression in scenarios.items():
            if scenario not in SCENARIOS:
                SCENARIOS.append(scenario)
            for sample in SAMPLES:
                for filetype in data["filetypes"]:
                    VALID_WILDCARD_COMBINATIONS.append(
                        {
                            "tool": tool,
                            "annotation": annotation,
                            "mode": "direct",
                            "scenario": scenario,
                            "sample": sample,
                            "filetype": (
                                filetype,
                                filetype[:3] + (".gz" if filetype == "vcf_z" else ""),
                            ),
                            "filecode": filetype,
                            "fileextension": filetype[:3]
                            + (".gz" if filetype == "vcf_z" else ""),
                        }
                    )

VALID_WILDCARD_COMBINATIONS = pd.DataFrame(VALID_WILDCARD_COMBINATIONS)


## Benchmark related constants and functions

BENCHMARK_REPEATS = int(config["benchmark"].get("repeats", 10))


def escape_expression(expression):
    return expression


## Validation related constants and functions


def group_validation_results(wildcards, *args, **kwargs):
    wc = wildcards
    tools = VALID_WILDCARD_COMBINATIONS.query(
        f"annotation == '{wc.annotation}' "
        f"and scenario == '{wc.scenario}' "
        f"and filecode == '{wc.filetype}'",
        engine="python",
    )["tool"].unique()
    paths = expand(
        "results/filtered/{annotation}/{tool}/{scenario}/{filetype}/{file}.sorted.noheader.md5sum",
        **wildcards,
        tool=tools,
    )
    return paths
