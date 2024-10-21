# statements
import pandas as pd
from pathlib import Path

def get_summary_bed(wildcards):
    for file in config["multiBamSummary"]["summary_bed"]:
        if Path(file).name.removesuffix(".bed") == wildcards.bed:
            return Path(file).absolute()

def get_rmsk():
    return Path("resources") / Path(config["rmsk_url"]).name

def get_genome():
    return Path("resources") / Path(config["genome_url"]).name

def get_effectiveGenomeSize(wildcards):
    if config["rmsk_url"]:
        df = pd.read_csv(Path("config") / "effectiveGenomeSize_rmsk.csv", header=0)
        return df.iloc[abs(df["Read length"] - config["fragmentLength"]).argsort()[0]][config["officialAssemblyName"]]
    df = pd.read_csv(Path("config") / "effectiveGenomeSize_normsk.csv", header=0)
    return df[df["Genome"] == config["officialAssemblyName"]]["Effective size"].values[0]

def get_sampleSize(wildcards, limit=None):
    return min(get_effectiveGenomeSize(wildcards), limit)

def get_sorted_bam(wildcards, suffix=".sorted.bam"):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    for file in df["sample"]:
        if Path(file).name.removesuffix(".sorted.bam") == wildcards.sample:
            return Path(file).absolute().as_posix().removesuffix(".sorted.bam") + suffix

def get_group_sorted_bam(wildcards, suffix=".sorted.bam"):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    return [
        Path("results") / (Path(file).name.removesuffix(".sorted.bam") + suffix)
        for file in df["sample"][df["group"] == wildcards.group]
    ]

def get_all_sorted_bam(wildcards, suffix=".sorted.bam", is_results=False):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    return [
        Path("results") / (Path(file).name.removesuffix(".sorted.bam") + suffix) if is_results
        else Path(file).absolute().as_posix().removesuffix(".sorted.bam") + suffix
        for file in df["sample"]
    ]

def get_all_samples(wildcards):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    return [
        Path(file).name.removesuffix(".sorted.bam")
        for file in df["sample"]
    ]

def get_all_groups(wildcards):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    return df["group"].unique()

def get_all_group_pairs(wildcards):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    group_pairs = list()
    for group in df["group"].unique():
        control_group = df["control_group"][df["group"] == group].unique()
        assert len(control_group) == 1, f"control_group for group {group} is not unique"
        control_group = control_group[0]
        if control_group:
            group_pairs.append([group, control_group])
    return group_pairs
        
def get_all_summary_bed_stem(wildcards):
    return [
        Path(file).stem
        for file in config["multiBamSummary"]["summary_bed"]
    ]

def get_all_bigwig(wildcards):
    samples = get_all_samples(wildcards)
    groups = get_all_groups(wildcards)
    group_pairs = get_all_group_pairs(wildcards)
    return [
        f'''results/{sample}.bigwig''' for sample in samples
    ] + [
        f'''results/{group}.averaged.bigwig''' for group in groups
    ] + [
        f'''results/{group}_cmp_{control_group}.bigwig''' for group, control_group in group_pairs
    ]

def get_all_bigwig_stem(wildcards):
    return [
        Path(file).stem for file in get_all_bigwig(wildcards)
    ]

def get_all_regionsFileName_stem(wildcards):
    return [
        Path(file).stem for file in config["computeMatrix"]["regionsFileName"]
    ]

def get_all_peaks(wildcards):
    df = pd.read_csv(Path("config") / "samplesheet.csv", header=0, na_filter=False)
    return df["peak"]

def get_run_all(wildcards):
    summary_bed_stems = get_all_summary_bed_stem(wildcards)
    return [
        "results/plotFingerprint.pdf",
        "results/bamPEFragmentSize.pdf",
        "results/plotCoverage.pdf",
        "results/plotEnrichment.pdf",
        "results/estimateReadFiltering.tsv",
        "results/computeMatrix.gz",
        "results/plotHeatmap.pdf"
    ] + [
        f"results/plotCorrelation.{stem}.pdf"
        for stem in summary_bed_stems
    ] + [
        f"results/plotPCA.{stem}.pdf"
        for stem in summary_bed_stems
    ]
