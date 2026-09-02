"""Build the exact submission layouts for Fig 1, Fig 6, and S6 Fig.

Scripts 07-16 perform the scientific analyses and create analytical panels.
This script is deliberately separate: it fixes the final panel geometry,
typography, labels, and TIFF export settings used for the PLOS ONE package.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


RESPONDER = "#F4A261"
NONRESPONDER = "#3DB7B3"
FIG1_SIZE = (4252, 3570)
FIG1_PANELS = {
    "Fig1A_layout.png": (0, 0),
    "Fig1B_layout.png": (2100, 0),
    "Fig1C_layout.png": (0, 1750),
    "Fig1D_layout.png": (2100, 1750),
}


def repository_root() -> Path:
    configured = os.environ.get("HCC_PROJECT_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(__file__).resolve().parents[1]


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        average = (start + 1 + stop) / 2.0
        for position in range(start, stop):
            ranks[order[position]] = average
        start = stop
    return ranks


def exact_rank_sum_p(x: list[float], y: list[float]) -> float:
    """Two-sided exact Wilcoxon rank-sum p value using all label permutations."""
    values = [float(value) for value in x] + [float(value) for value in y]
    n_x = len(x)
    ranks = average_ranks(values)
    observed = sum(ranks[:n_x])
    permutations = [
        sum(ranks[index] for index in combination)
        for combination in itertools.combinations(range(len(values)), n_x)
    ]
    lower = sum(value <= observed + 1e-12 for value in permutations) / len(permutations)
    upper = sum(value >= observed - 1e-12 for value in permutations) / len(permutations)
    return min(1.0, 2.0 * min(lower, upper))


def p_label(value: float, digits: int = 3) -> str:
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.{digits}f}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_tiff(image: Image.Image, path: Path, dpi: int = 600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").save(path, compression="tiff_lzw", dpi=(dpi, dpi))


def save_matplotlib_figure(fig, tif_path: Path, pdf_path: Path, dpi: int = 600) -> None:
    tif_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, facecolor="white")
    preview = tif_path.with_suffix(".build-preview.png")
    fig.savefig(preview, dpi=dpi, facecolor="white")
    with Image.open(preview) as image:
        save_tiff(image, tif_path, dpi=dpi)
    preview.unlink()


def build_fig1(assets: Path, output: Path) -> Path:
    """Compose the locked submission geometry from versioned panel exports."""
    canvas = Image.new("RGB", FIG1_SIZE, "white")
    for filename, position in FIG1_PANELS.items():
        panel_path = assets / filename
        if not panel_path.exists():
            raise FileNotFoundError(f"Missing Fig 1 layout asset: {panel_path}")
        with Image.open(panel_path) as panel:
            canvas.paste(panel.convert("RGB"), position)
    destination = output / "Fig1.tif"
    save_tiff(canvas, destination)
    return destination


def boxplot_panel(
    ax,
    responder: list[float],
    nonresponder: list[float],
    title: str,
    subtitle: str,
    ylabel: str,
    p_value: float,
    panel: str,
    digits: int = 3,
    xlabels: tuple[str, str] = ("Responder", "Nonresponder"),
) -> None:
    groups = [responder, nonresponder]
    bp = ax.boxplot(
        groups,
        positions=[1, 2],
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#2F2F2F", "linewidth": 1.3},
        whiskerprops={"color": "#404040", "linewidth": 0.9},
        capprops={"color": "#404040", "linewidth": 0.9},
        boxprops={"edgecolor": "#404040", "linewidth": 0.9},
    )
    for patch, color in zip(bp["boxes"], (RESPONDER, NONRESPONDER)):
        patch.set_facecolor(color)
        patch.set_alpha(0.88)
    for index, values in enumerate(groups, start=1):
        offsets = np.linspace(-0.075, 0.075, len(values)) if len(values) > 1 else np.array([0.0])
        ax.scatter(index + offsets, values, s=19, color="#222222", zorder=3, edgecolor="none")
    combined = np.asarray(responder + nonresponder, dtype=float)
    span = max(float(np.ptp(combined)), 0.05)
    ax.set_ylim(float(combined.min() - 0.18 * span), float(combined.max() + 0.33 * span))
    ax.text(1.5, float(combined.max() + 0.20 * span), p_label(p_value, digits), ha="center", fontsize=8)
    ax.set_xticks([1, 2], xlabels)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_title(title, fontsize=9.3, pad=28)
    ax.text(
        0.5,
        1.015,
        subtitle,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
        color="#777777",
        style="italic",
        linespacing=1.05,
    )
    ax.text(-0.12, 1.10, panel, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def forest_panel(ax, data: pd.DataFrame, variable: str, title: str, subtitle: str, color: str, panel: str) -> None:
    cohort_order = ["GSE140901", "GSE279750", "GSE215011"]
    subset = data[data["variable"] == variable].copy().set_index("dataset").loc[cohort_order].reset_index()
    y = np.array([2, 1, 0], dtype=float)
    estimate = subset["cohens_d"].to_numpy(float)
    low = subset["d_CI_lower"].to_numpy(float)
    high = subset["d_CI_upper"].to_numpy(float)
    ax.axvline(0, color="#9A9A9A", linestyle="--", linewidth=0.9, zorder=0)
    ax.errorbar(
        estimate,
        y,
        xerr=np.vstack([estimate - low, high - estimate]),
        fmt="o",
        color=color,
        ecolor="#303030",
        elinewidth=1.1,
        capsize=5,
        markersize=6.2,
        zorder=3,
    )
    treatments = {
        "GSE140901": "anti-PD-1",
        "GSE279750": "anti-PD-L1",
        "GSE215011": "anti-PD-1",
    }
    ax.set_yticks(y, [f"{cohort}\n{treatments[cohort]}" for cohort in subset["dataset"]])
    for ypos, effect, p_value in zip(y, estimate, subset["wilcox_p"].to_numpy(float)):
        ax.text(
            effect + 0.16,
            ypos + 0.12,
            f"d={effect:.2f}; {p_label(p_value)}",
            va="bottom",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6, "alpha": 0.9},
        )
    ax.set_xlim(min(-2.7, float(low.min() - 0.25)), max(3.2, float(high.max() + 0.55)))
    ax.set_ylim(-0.55, 2.55)
    ax.set_xlabel("Cohen's d (responder vs nonresponder)", fontsize=8.5)
    ax.set_title(title, fontsize=9.3, pad=28)
    ax.text(0.5, 1.015, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=8, color="#777777", style="italic")
    ax.text(-0.12, 1.10, panel, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def build_fig6(data_dir: Path, output: Path) -> tuple[Path, dict[str, float]]:
    gse140 = pd.read_csv(data_dir / "GSE140901_scores_normalized_method.csv")
    gse140 = gse140[gse140["response"].isin(["Responder", "NonResponder"])]
    gse279 = pd.read_csv(data_dir / "GSE279750_signature_scores.csv")
    cross = pd.read_csv(data_dir / "cross_cohort_effect_sizes.csv")

    gse140_r = gse140.loc[gse140["response"] == "Responder", "cxcl13_score"].tolist()
    gse140_nr = gse140.loc[gse140["response"] == "NonResponder", "cxcl13_score"].tolist()
    gse279_r = gse279.loc[gse279["response"] == "Responder", "cxcl13_associated_exhaustion_score"].tolist()
    gse279_nr = gse279.loc[gse279["response"] == "Nonresponder", "cxcl13_associated_exhaustion_score"].tolist()
    p_140 = exact_rank_sum_p(gse140_r, gse140_nr)
    p_279 = exact_rank_sum_p(gse279_r, gse279_nr)

    fig = plt.figure(figsize=(7.5, 6.5), constrained_layout=False)
    grid = fig.add_gridspec(2, 2, left=0.12, right=0.98, bottom=0.09, top=0.86, hspace=0.70, wspace=0.62)
    boxplot_panel(
        fig.add_subplot(grid[0, 0]),
        gse140_r,
        gse140_nr,
        "GSE140901: CXCL13-associated\nexhaustion score",
        "PR=6, PD=8 | clinical RECIST labels\nNanoString-normalization sensitivity analysis",
        "CXCL13-associated exhaustion score",
        p_140,
        "A",
    )
    boxplot_panel(
        fig.add_subplot(grid[0, 1]),
        gse279_r,
        gse279_nr,
        "GSE279750: CXCL13-associated\nexhaustion score",
        "anti-PD-L1 | 6 responders vs 4 nonresponders\nbiopsy timepoint not explicitly stated",
        "CXCL13-associated exhaustion score",
        p_279,
        "B",
    )
    forest_panel(
        fig.add_subplot(grid[1, 0]),
        cross,
        "cxcl13_associated_exhaustion_score",
        "CXCL13-associated exhaustion score\nacross clinical ICI cohorts",
        "Small cohorts; effect directions are exploratory",
        "#E41A1C",
        "C",
    )
    forest_panel(
        fig.add_subplot(grid[1, 1]),
        cross,
        "tpex_score",
        "Tpex-like score across\nclinical ICI cohorts",
        "Direction is inconsistent across cohorts",
        "#377EB8",
        "D",
    )
    destination = output / "Fig6.tif"
    save_matplotlib_figure(fig, destination, output / "Fig6.pdf")
    plt.close(fig)
    return destination, {"GSE140901_exact_p": p_140, "GSE279750_exact_p": p_279}


def cox_panel(ax, data: pd.DataFrame) -> None:
    labels = []
    for _, row in data.iterrows():
        score = {"tpex_score": "Tpex-like", "cxcl13_score": "CXCL13", "combined_score": "Combined"}[row["variable"]]
        model_id = "M2" if "grade" in row["model"] else "M1"
        labels.append(f"{score} ({model_id})")
    y = np.arange(len(data))[::-1]
    hr = data["HR"].to_numpy(float)
    low = data["HR_lower"].to_numpy(float)
    high = data["HR_upper"].to_numpy(float)
    ax.axvline(1, color="#999999", linestyle="--", linewidth=0.9)
    ax.errorbar(hr, y, xerr=np.vstack([hr - low, high - hr]), fmt="o", color="#9C9C9C", ecolor="#333333", capsize=4, markersize=5)
    for ypos, estimate, p_value in zip(y, hr, data["p"].to_numpy(float)):
        ax.text(estimate + 0.08, ypos + 0.11, p_label(p_value), va="bottom", fontsize=8, bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4, "alpha": 0.9})
    ax.set_yticks(y, labels)
    ax.set_xlim(-0.15, 3.25)
    ax.set_xlabel("Hazard ratio (95% CI)", fontsize=8)
    ax.set_title("Multivariable Cox sensitivity\nM1: age/sex/stage; M2: + grade", fontsize=9, pad=14, linespacing=1.15)
    ax.text(-0.18, 1.16, "A", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def normalization_panel(ax, log2_scores: pd.DataFrame, norm_scores: pd.DataFrame) -> tuple[float, float]:
    log2_scores = log2_scores[log2_scores["response"].isin(["Responder", "NonResponder"])]
    norm_scores = norm_scores[norm_scores["response"].isin(["Responder", "NonResponder"])]
    groups = [
        log2_scores.loc[log2_scores["response"] == "Responder", "cxcl13_score"].tolist(),
        log2_scores.loc[log2_scores["response"] == "NonResponder", "cxcl13_score"].tolist(),
        norm_scores.loc[norm_scores["response"] == "Responder", "cxcl13_score"].tolist(),
        norm_scores.loc[norm_scores["response"] == "NonResponder", "cxcl13_score"].tolist(),
    ]
    p_log2 = exact_rank_sum_p(groups[0], groups[1])
    p_norm = exact_rank_sum_p(groups[2], groups[3])
    positions = [1, 2, 4, 5]
    bp = ax.boxplot(
        groups,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#333333", "linewidth": 1.1},
        whiskerprops={"color": "#444444", "linewidth": 0.8},
        capprops={"color": "#444444", "linewidth": 0.8},
        boxprops={"edgecolor": "#444444", "linewidth": 0.8},
    )
    for patch, color in zip(bp["boxes"], [RESPONDER, NONRESPONDER, RESPONDER, NONRESPONDER]):
        patch.set_facecolor(color)
        patch.set_alpha(0.88)
    for position, values in zip(positions, groups):
        offsets = np.linspace(-0.06, 0.06, len(values))
        ax.scatter(position + offsets, values, s=12, color="#222222", zorder=3)
    ax.axvline(3, color="#CCCCCC", linewidth=0.8)
    ax.set_xticks(positions, ["Responder\nlog2", "Nonresp.\nlog2", "Responder\nNanoString", "Nonresp.\nNanoString"])
    ymax = max(max(group) for group in groups)
    ymin = min(min(group) for group in groups)
    span = ymax - ymin
    ax.set_ylim(ymin - 0.12 * span, ymax + 0.30 * span)
    ax.text(1.5, ymax + 0.16 * span, p_label(p_log2), ha="center", fontsize=8)
    ax.text(4.5, ymax + 0.16 * span, p_label(p_norm), ha="center", fontsize=8)
    ax.set_ylabel("CXCL13-associated score", fontsize=8)
    ax.set_title("GSE140901 normalization sensitivity", fontsize=9, pad=14)
    ax.text(-0.12, 1.16, "B", transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    return p_log2, p_norm


def build_s6(data_dir: Path, output: Path) -> tuple[Path, dict[str, float]]:
    cox = pd.read_csv(data_dir / "TCGA_multivariate_sensitivity.csv")
    log2_scores = pd.read_csv(data_dir / "GSE140901_scores_original_method.csv")
    norm_scores = pd.read_csv(data_dir / "GSE140901_scores_normalized_method.csv")
    fad = pd.read_csv(data_dir / "GSE202069_rerun_scores.csv")
    gse215 = pd.read_csv(data_dir / "GSE215011_signature_scores.csv")

    fig = plt.figure(figsize=(7.5, 8.5), constrained_layout=False)
    grid = fig.add_gridspec(3, 12, left=0.13, right=0.985, bottom=0.06, top=0.90, hspace=0.82, wspace=1.05)
    cox_panel(fig.add_subplot(grid[0, 0:6]), cox)
    p_log2, p_norm = normalization_panel(fig.add_subplot(grid[0, 6:12]), log2_scores, norm_scores)

    panels = [
        (grid[1, 0:4], "tpex_score", "Tpex-like score", "Tpex-like score", "C", 4),
        (grid[1, 4:8], "cxcl13_score", "CXCL13-associated\nexhaustion score", "CXCL13-associated score", "D", 3),
        (grid[1, 8:12], "combined_score", "Combined score", "Combined score", "E", 2),
    ]
    statistics: dict[str, float] = {"GSE140901_log2_exact_p": p_log2, "GSE140901_NanoString_exact_p": p_norm}
    for spec, column, title, ylabel, letter, digits in panels:
        responder = fad.loc[fad["group"] == "FAD_Responder", column].tolist()
        nonresponder = fad.loc[fad["group"] == "FAD_NonResponder", column].tolist()
        p_value = exact_rank_sum_p(responder, nonresponder)
        statistics[f"GSE202069_{column}_exact_p"] = p_value
        boxplot_panel(
            fig.add_subplot(spec),
            responder,
            nonresponder,
            title,
            "GSE202069 | FAD-inferred\nsubgroups (n=6 vs 5)",
            ylabel,
            p_value,
            letter,
            digits,
            ("FAD-inferred R", "FAD-inferred NR"),
        )

    external = [
        (grid[2, 1:6], "cxcl13_associated_exhaustion_score", "CXCL13-associated\nexhaustion score", "CXCL13-associated score", "F", 2),
        (grid[2, 6:11], "tpex_score", "Tpex-like score", "Tpex-like score", "G", 2),
    ]
    for spec, column, title, ylabel, letter, digits in external:
        responder = gse215.loc[gse215["response"] == "Responder", column].tolist()
        nonresponder = gse215.loc[gse215["response"] == "Nonresponder", column].tolist()
        p_value = exact_rank_sum_p(responder, nonresponder)
        statistics[f"GSE215011_{column}_exact_p"] = p_value
        boxplot_panel(
            fig.add_subplot(spec),
            responder,
            nonresponder,
            title,
            "GSE215011 | nivolumab\n5 responders vs 5 nonresponders",
            ylabel,
            p_value,
            letter,
            digits,
        )

    destination = output / "S6_Fig.tif"
    save_matplotlib_figure(fig, destination, output / "S6_Fig.pdf")
    plt.close(fig)
    return destination, statistics


def inspect_tiff(path: Path) -> dict[str, object]:
    with Image.open(path) as image:
        dpi = image.info.get("dpi", (None, None))
        return {
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "pixels": list(image.size),
            "dpi": [float(value) if value is not None else None for value in dpi],
            "mode": image.mode,
            "compression": image.info.get("compression"),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Override the default submission_figures directory")
    args = parser.parse_args()

    root = repository_root()
    data_dir = root / "processed_data"
    assets = root / "figure_assets"
    output = args.output_dir.resolve() if args.output_dir else root / "submission_figures"
    output.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": os.environ.get("HCC_FIGURE_FONT", "Arial"),
            "font.size": 8,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig1 = build_fig1(assets, output)
    fig6, fig6_stats = build_fig6(data_dir, output)
    s6, s6_stats = build_s6(data_dir, output)
    manifest = {
        "purpose": "Versioned construction of the final PLOS ONE submission layouts",
        "repository_root": ".",
        "font_family": plt.rcParams["font.family"],
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "figures": {path.name: inspect_tiff(path) for path in (fig1, fig6, s6)},
        "fig1_layout_assets": {name: sha256(assets / name) for name in FIG1_PANELS},
        "displayed_statistics": {**fig6_stats, **s6_stats},
    }
    manifest_path = output / "layout_build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
