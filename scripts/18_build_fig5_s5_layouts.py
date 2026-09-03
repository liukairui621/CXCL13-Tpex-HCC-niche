"""Build the final PLOS ONE layouts for Fig 5 and S5 Fig.

The script reads only versioned processed tables, preserves the submitted
statistics, and controls panel geometry, typography, and 600-dpi TIFF export.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import platform
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from PIL import Image
import scipy
from scipy import stats


RED = "#E41A1C"
BLUE = "#377EB8"
GREEN = "#4DAF4A"
BLACK = "#222222"

FONT_SIZES = {
    "tick": 8.0,
    "axis": 8.7,
    "legend": 8.0,
    "annotation": 8.0,
    "title": 9.2,
    "panel": 12.0,
}


def repository_root() -> Path:
    configured = os.environ.get("HCC_PROJECT_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": os.environ.get("HCC_FIGURE_FONT", "Arial"),
            "font.size": FONT_SIZES["tick"],
            "axes.labelsize": FONT_SIZES["axis"],
            "axes.titlesize": FONT_SIZES["title"],
            "xtick.labelsize": FONT_SIZES["tick"],
            "ytick.labelsize": FONT_SIZES["tick"],
            "legend.fontsize": FONT_SIZES["legend"],
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def style_axis(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(BLACK)
    ax.spines["bottom"].set_color(BLACK)
    ax.tick_params(axis="both", colors=BLACK, labelsize=FONT_SIZES["tick"], length=3)


def save_figure(fig: mpl.figure.Figure, tif_path: Path, pdf_path: Path, dpi: int = 600) -> None:
    tif_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, bbox_inches=None, facecolor="white")
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches=None, facecolor="white")
    buffer.seek(0)
    with Image.open(buffer) as source:
        source.convert("RGB").save(
            tif_path,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(dpi, dpi),
        )
    plt.close(fig)


def km_estimate(durations: np.ndarray, events: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    durations = np.asarray(durations, dtype=float)
    events = np.asarray(events, dtype=int)
    event_times = np.sort(np.unique(durations[events == 1]))
    times = [0.0]
    survival = [1.0]
    current = 1.0
    for time in event_times:
        at_risk = np.sum(durations >= time)
        deaths = np.sum((durations == time) & (events == 1))
        if at_risk:
            current *= 1.0 - deaths / at_risk
            times.append(float(time))
            survival.append(float(current))
    max_time = float(np.nanmax(durations))
    if times[-1] < max_time:
        times.append(max_time)
        survival.append(current)
    return np.asarray(times), np.asarray(survival)


def km_survival_at(durations: np.ndarray, events: np.ndarray, query: np.ndarray) -> np.ndarray:
    times, survival = km_estimate(durations, events)
    indices = np.searchsorted(times, query, side="right") - 1
    indices = np.clip(indices, 0, len(survival) - 1)
    return survival[indices]


def multigroup_logrank(durations: np.ndarray, events: np.ndarray, groups: np.ndarray) -> tuple[float, float]:
    group_levels = list(pd.unique(groups))
    k = len(group_levels)
    event_times = np.sort(np.unique(durations[events == 1]))
    observed = np.zeros(k)
    expected = np.zeros(k)
    covariance = np.zeros((k, k))

    for time in event_times:
        risk = np.asarray([np.sum((groups == level) & (durations >= time)) for level in group_levels], dtype=float)
        deaths = np.asarray(
            [np.sum((groups == level) & (durations == time) & (events == 1)) for level in group_levels],
            dtype=float,
        )
        total_risk = risk.sum()
        total_deaths = deaths.sum()
        if total_risk <= 1 or total_deaths == 0:
            continue
        observed += deaths
        expected += total_deaths * risk / total_risk
        factor = total_deaths * (total_risk - total_deaths) / (total_risk - 1)
        proportions = risk / total_risk
        covariance += factor * (np.diag(proportions) - np.outer(proportions, proportions))

    delta = (observed - expected)[:-1]
    covariance_reduced = covariance[:-1, :-1]
    statistic = float(delta.T @ np.linalg.pinv(covariance_reduced) @ delta)
    p_value = float(stats.chi2.sf(statistic, k - 1))
    return statistic, p_value


def p_text(value: float) -> str:
    if value < 0.001:
        return "p < 0.001"
    if value < 0.01:
        return f"p = {value:.4f}"
    return f"p = {value:.2f}"


def add_panel_letter(ax: mpl.axes.Axes, letter: str, x: float = -0.18, y: float = 1.15) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=FONT_SIZES["panel"],
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def draw_km_panel(
    main_ax: mpl.axes.Axes,
    risk_title_ax: mpl.axes.Axes,
    risk_ax: mpl.axes.Axes,
    data: pd.DataFrame,
    group_col: str,
    group_order: list[str],
    colors: dict[str, str],
    title: str,
    panel_letter: str,
) -> float:
    durations = data["OS.time"].to_numpy(float)
    events = data["OS"].to_numpy(int)
    groups = data[group_col].to_numpy(str)
    _, p_value = multigroup_logrank(durations, events, groups)
    max_display = 4200

    for group in group_order:
        subset = data[data[group_col] == group]
        group_durations = subset["OS.time"].to_numpy(float)
        group_events = subset["OS"].to_numpy(int)
        times, survival = km_estimate(group_durations, group_events)
        main_ax.step(times, survival, where="post", color=colors[group], linewidth=1.5, label=group)
        censored = group_durations[group_events == 0]
        if len(censored):
            censored_survival = km_survival_at(group_durations, group_events, censored)
            main_ax.plot(
                censored,
                censored_survival,
                linestyle="none",
                marker="+",
                markersize=4.0,
                markeredgewidth=0.75,
                color=colors[group],
            )

    style_axis(main_ax)
    main_ax.set_xlim(-100, max_display)
    main_ax.set_ylim(-0.02, 1.05)
    main_ax.set_xticks([0, 1000, 2000, 3000, 4000])
    main_ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    main_ax.set_xlabel("")
    main_ax.set_ylabel("Overall survival probability", labelpad=4)
    main_ax.set_title(title, y=1.17, pad=0, fontsize=FONT_SIZES["title"])
    main_ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.105),
        ncol=len(group_order),
        frameon=False,
        handlelength=1.2,
        columnspacing=0.9,
        handletextpad=0.35,
    )
    main_ax.text(
        0.055,
        0.18,
        p_text(p_value),
        transform=main_ax.transAxes,
        fontsize=FONT_SIZES["annotation"],
    )
    add_panel_letter(main_ax, panel_letter)

    risk_title_ax.axis("off")
    risk_title_ax.text(
        0.5,
        0.25,
        "Time (days)",
        transform=risk_title_ax.transAxes,
        ha="center",
        va="center",
        fontsize=FONT_SIZES["axis"],
    )
    risk_ax.set_xlim(-100, max_display)
    risk_ax.set_ylim(0, 1)
    risk_ax.axis("off")
    risk_ax.text(
        0.5,
        0.82,
        "Number at risk",
        transform=risk_ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=FONT_SIZES["axis"],
        clip_on=False,
    )
    risk_times = [0, 1000, 2000, 3000, 4000]
    if len(group_order) == 2:
        row_positions = [0.45, 0.08]
    else:
        row_positions = [0.55, 0.30, 0.05]
    for y, group in zip(row_positions, group_order):
        subset = data[data[group_col] == group]
        risk_ax.text(
            -0.035,
            y,
            group,
            transform=risk_ax.transAxes,
            clip_on=False,
            color=colors[group],
            ha="right",
            va="center",
            fontsize=FONT_SIZES["legend"],
        )
        for time in risk_times:
            count = int(np.sum(subset["OS.time"].to_numpy(float) >= time))
            risk_ax.text(time, y, str(count), ha="center", va="center", fontsize=FONT_SIZES["tick"], color=BLACK)

    return p_value


def build_fig5(data_dir: Path, output_dir: Path) -> dict[str, float]:
    source = pd.read_csv(data_dir / "TCGA_LIHC_scores_and_survival.csv")
    # Preserve the submitted analysis order: define median-split groups in the
    # complete score table, then exclude the one record lacking OS follow-up.
    source["tpex_group"] = np.where(source["tpex_score"] > source["tpex_score"].median(), "Tpex High", "Tpex Low")
    source["cxcl13_group"] = np.where(
        source["cxcl13_associated_exhaustion_score"] > source["cxcl13_associated_exhaustion_score"].median(),
        "CXCL13 High",
        "CXCL13 Low",
    )
    source["combined_group"] = np.where(
        (source["tpex_group"] == "Tpex High") & (source["cxcl13_group"] == "CXCL13 High"),
        "Both High",
        np.where(
            (source["tpex_group"] == "Tpex Low") & (source["cxcl13_group"] == "CXCL13 Low"),
            "Both Low",
            "Mixed",
        ),
    )
    source = source.dropna(subset=["OS.time", "OS", "tpex_score", "cxcl13_associated_exhaustion_score"]).copy()
    source["OS"] = source["OS"].astype(int)

    fig = plt.figure(figsize=(7.5, 8.25))
    # Keep the four-panel structure while reserving enough horizontal space
    # for the forest-plot row labels so they cannot intrude into panel C.
    outer = fig.add_gridspec(
        2,
        2,
        left=0.105,
        right=0.985,
        top=0.905,
        bottom=0.055,
        wspace=0.52,
        hspace=0.38,
    )
    panel_specs = [
        (outer[0, 0], "tpex_group", ["Tpex High", "Tpex Low"], {"Tpex High": RED, "Tpex Low": BLUE}, "Tpex-like signature (TCGA-LIHC)", "A"),
        (outer[0, 1], "cxcl13_group", ["CXCL13 High", "CXCL13 Low"], {"CXCL13 High": RED, "CXCL13 Low": BLUE}, "CXCL13-associated exhaustion (TCGA-LIHC)", "B"),
        (outer[1, 0], "combined_group", ["Both High", "Mixed", "Both Low"], {"Both High": RED, "Mixed": GREEN, "Both Low": BLUE}, "Combined signature (TCGA-LIHC)", "C"),
    ]

    results: dict[str, float] = {}
    for spec, group_col, group_order, colors, title, letter in panel_specs:
        nested = spec.subgridspec(3, 1, height_ratios=[2.75, 0.65, 0.95], hspace=0.00)
        main_ax = fig.add_subplot(nested[0])
        risk_title_ax = fig.add_subplot(nested[1])
        risk_ax = fig.add_subplot(nested[2])
        results[f"Fig5{letter}_logrank_p"] = draw_km_panel(
            main_ax, risk_title_ax, risk_ax, source, group_col, group_order, colors, title, letter
        )

    d_nested = outer[1, 1].subgridspec(3, 1, height_ratios=[2.75, 0.65, 0.95], hspace=0.00)
    ax_d = fig.add_subplot(d_nested[0])
    fig.add_subplot(d_nested[1]).axis("off")
    fig.add_subplot(d_nested[2]).axis("off")
    cox = pd.read_csv(data_dir / "TCGA_cox_results_continuous.csv")
    label_map = {
        "tpex_score": "Tpex-like score",
        "cxcl13_score": "CXCL13-associated\nscore",
        "combined_score": "Combined score",
    }
    order = ["tpex_score", "cxcl13_score", "combined_score"]
    cox = cox.set_index("variable").loc[order].reset_index()
    y_positions = np.arange(len(cox))[::-1]
    ax_d.axvline(1.0, color="#999999", linestyle="--", linewidth=0.8)
    for y, row in zip(y_positions, cox.itertuples(index=False)):
        color = RED if row.p < 0.05 else "#999999"
        ax_d.errorbar(
            row.HR,
            y,
            xerr=[[row.HR - row.HR_lower], [row.HR_upper - row.HR]],
            fmt="o",
            markersize=5,
            color=color,
            ecolor=BLACK,
            elinewidth=1.1,
            capsize=4,
            capthick=1.1,
            zorder=3,
        )
        ax_d.text(
            row.HR + 0.055,
            y + 0.17,
            f"HR={row.HR:.2f}, p={row.p:.3f}",
            fontsize=FONT_SIZES["annotation"],
            ha="left",
            va="bottom",
            color=BLACK,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.35, "alpha": 0.92},
        )
    style_axis(ax_d)
    ax_d.set_xlim(-0.05, 2.55)
    ax_d.set_ylim(-0.55, 2.65)
    ax_d.set_yticks(y_positions)
    ax_d.set_yticklabels([label_map[value] for value in order])
    ax_d.set_xlabel("Hazard ratio (95% CI)")
    ax_d.set_title("Univariate Cox regression (TCGA-LIHC)", pad=28)
    ax_d.text(
        0.5,
        1.03,
        f"Continuous scores, N={int(cox['n'].iloc[0])}, events={int(cox['n_events'].iloc[0])}",
        transform=ax_d.transAxes,
        ha="center",
        va="bottom",
        fontsize=FONT_SIZES["annotation"],
        color="#777777",
        style="italic",
    )
    add_panel_letter(ax_d, "D", x=-0.15, y=1.15)

    fig5_path = output_dir / "Fig5.tif"
    save_figure(fig, fig5_path, output_dir / "Fig5.pdf")
    results["Fig5_n"] = float(len(source))
    results["Fig5_file"] = inspect_tiff(fig5_path)
    return results


def add_linear_fit(ax: mpl.axes.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    fit = stats.linregress(x, y)
    x_grid = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), 160)
    y_grid = fit.intercept + fit.slope * x_grid
    residuals = y - (fit.intercept + fit.slope * x)
    degrees = max(len(x) - 2, 1)
    residual_se = math.sqrt(float(np.sum(residuals**2)) / degrees)
    sxx = float(np.sum((x - np.mean(x)) ** 2))
    if sxx > 0:
        mean_se = residual_se * np.sqrt(1.0 / len(x) + (x_grid - np.mean(x)) ** 2 / sxx)
        critical = stats.t.ppf(0.975, degrees)
        ax.fill_between(x_grid, y_grid - critical * mean_se, y_grid + critical * mean_se, color=color, alpha=0.13, linewidth=0)
    ax.plot(x_grid, y_grid, color=color, linewidth=1.4)


def format_rho(value: float) -> str:
    return f"{value:.2f}" if abs(value) >= 0.1 else f"{value:.3f}"


def draw_spot_violin(
    ax: mpl.axes.Axes,
    data: pd.DataFrame,
    y_col: str,
    y_label: str,
    title: str,
    letter: str,
    random_seed: int,
) -> float:
    order = ["Responder", "NonResponder"]
    palette = {"Responder": RED, "NonResponder": BLUE}
    sns.violinplot(
        data=data,
        x="response",
        y=y_col,
        order=order,
        hue="response",
        palette=palette,
        legend=False,
        inner=None,
        cut=0,
        linewidth=0.7,
        density_norm="width",
        ax=ax,
    )
    sns.boxplot(
        data=data,
        x="response",
        y=y_col,
        order=order,
        width=0.22,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": BLACK, "linewidth": 0.7, "zorder": 4},
        whiskerprops={"color": BLACK, "linewidth": 0.7},
        capprops={"color": BLACK, "linewidth": 0.7},
        medianprops={"color": BLACK, "linewidth": 0.9},
        ax=ax,
    )
    rng = np.random.default_rng(random_seed)
    for position, group in enumerate(order):
        values = data.loc[data["response"] == group, y_col].to_numpy(float)
        jitter = rng.normal(position, 0.035, size=len(values))
        ax.scatter(jitter, values, s=1.2, color=BLACK, alpha=0.22, linewidths=0, rasterized=True, zorder=2)

    responder = data.loc[data["response"] == "Responder", y_col]
    nonresponder = data.loc[data["response"] == "NonResponder", y_col]
    p_value = float(stats.mannwhitneyu(responder, nonresponder, alternative="two-sided", method="asymptotic").pvalue)
    style_axis(ax)
    ax.set_xlabel("")
    ax.set_ylabel(y_label)
    ax.set_title(title, pad=9, linespacing=1.12)
    ax.text(
        0.5,
        0.94,
        "Descriptive Wilcoxon\n" + p_text(p_value),
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=FONT_SIZES["annotation"],
        linespacing=1.05,
    )
    # Keep the panel letter between the y-axis spine and the centered title;
    # long vertical labels otherwise collide with letters placed farther left.
    add_panel_letter(ax, letter, x=-0.13, y=1.15)
    return p_value


def draw_correlation(
    ax: mpl.axes.Axes,
    data: pd.DataFrame,
    y_col: str,
    y_label: str,
    title: str,
    letter: str,
    annotation_location: str,
) -> dict[str, float]:
    colors = {"Responder": RED, "NonResponder": BLUE}
    results: dict[str, float] = {}
    for group in ["Responder", "NonResponder"]:
        subset = data[data["response"] == group]
        x = subset["tpex_score"].to_numpy(float)
        y = subset[y_col].to_numpy(float)
        ax.scatter(x, y, s=2.1, color=colors[group], alpha=0.28, linewidths=0, rasterized=True)
        add_linear_fit(ax, x, y, colors[group])
        rho, p_value = stats.spearmanr(x, y)
        results[f"{group}_rho"] = float(rho)
        results[f"{group}_p"] = float(p_value)

    style_axis(ax)
    ax.set_xlabel("Tpex-like score")
    ax.set_ylabel(y_label)
    ax.set_title(title, pad=8)
    labels = []
    for group in ["Responder", "NonResponder"]:
        rho = results[f"{group}_rho"]
        p_value = results[f"{group}_p"]
        labels.append((colors[group], f"R = {format_rho(rho)}, {p_text(p_value)}"))
    if annotation_location == "top_left":
        positions = [(0.06, 0.94), (0.06, 0.86)]
        va = "top"
        ha = "left"
    else:
        positions = [(0.97, 0.17), (0.97, 0.09)]
        va = "bottom"
        ha = "right"
    for (color, label), position in zip(labels, positions):
        ax.text(
            *position,
            label,
            transform=ax.transAxes,
            color=color,
            fontsize=FONT_SIZES["annotation"],
            va=va,
            ha=ha,
            fontstyle="italic",
        )
    add_panel_letter(ax, letter, x=-0.18, y=1.14)
    return results


def build_s5(data_dir: Path, output_dir: Path) -> dict[str, float]:
    data = pd.read_csv(data_dir / "Q3_all_scores.csv")
    data["response"] = pd.Categorical(data["response"], categories=["Responder", "NonResponder"], ordered=True)
    data["combined_score"] = data["tpex_score"] + data["cxcl13_score"] + data["tls_score"]

    # Four violin plots are unreadably narrow in a single 7.5-inch row.
    # A 2 + 2 + 3 layout keeps the panel order but gives every y-axis title,
    # panel letter, and plot title an independent margin.
    fig = plt.figure(figsize=(7.5, 8.5))
    outer = fig.add_gridspec(
        3,
        1,
        left=0.10,
        right=0.985,
        top=0.925,
        bottom=0.105,
        hspace=0.76,
        height_ratios=[1.0, 1.0, 1.10],
    )
    top_grid = outer[0].subgridspec(1, 2, wspace=0.42)
    middle_grid = outer[1].subgridspec(1, 2, wspace=0.42)
    bottom_grid = outer[2].subgridspec(1, 3, wspace=0.58)
    top_axes = [fig.add_subplot(top_grid[0, index]) for index in range(2)]
    middle_axes = [fig.add_subplot(middle_grid[0, index]) for index in range(2)]
    violin_axes = top_axes + middle_axes
    bottom_axes = [fig.add_subplot(bottom_grid[0, index]) for index in range(3)]
    results: dict[str, float] = {}
    violin_specs = [
        ("tpex_score", "Tpex-like score", "Tpex-like signature", "A"),
        ("cxcl13_score", "CXCL13-associated exhaustion score", "CXCL13-associated\nexhaustion", "B"),
        ("tls_score", "TLS score", "TLS signature", "C"),
        ("combined_score", "Combined score", "Combined niche score", "D"),
    ]
    for index, (axis, spec) in enumerate(zip(violin_axes, violin_specs)):
        y_col, y_label, title, letter = spec
        results[f"S5{letter}_wilcoxon_p"] = draw_spot_violin(
            axis, data, y_col, y_label, title, letter, random_seed=20260902 + index
        )

    correlation_specs = [
        (
            "cxcl13_score",
            "CXCL13-associated exhaustion score",
            "Tpex-like vs CXCL13-associated\nexhaustion",
            "E",
            "top_left",
        ),
        ("tls_score", "TLS score", "Tpex-like vs TLS", "F", "top_left"),
        ("spp1_score", "SPP1 suppressive score", "Tpex-like vs SPP1 niche", "G", "bottom_right"),
    ]
    for axis, spec in zip(bottom_axes, correlation_specs):
        y_col, y_label, title, letter, location = spec
        panel_results = draw_correlation(axis, data, y_col, y_label, title, letter, location)
        for key, value in panel_results.items():
            results[f"S5{letter}_{key}"] = value

    legend_handles = [
        Line2D([0], [0], marker="o", color=color, markerfacecolor=color, markersize=4, linewidth=1.2, label=group)
        for group, color in [("Responder", RED), ("NonResponder", BLUE)]
    ]
    fig.legend(
        handles=legend_handles,
        title="response",
        title_fontsize=FONT_SIZES["legend"],
        loc="center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=2,
        frameon=False,
        handlelength=1.2,
        columnspacing=1.2,
        handletextpad=0.4,
    )

    s5_path = output_dir / "S5_Fig.tif"
    save_figure(fig, s5_path, output_dir / "S5_Fig.pdf")

    panel_fig, panel_ax = plt.subplots(figsize=(6, 5))
    panel_fig.subplots_adjust(left=0.16, right=0.97, top=0.83, bottom=0.14)
    draw_correlation(
        panel_ax,
        data,
        "cxcl13_score",
        "CXCL13-associated exhaustion score",
        "Tpex-like vs CXCL13-associated\nexhaustion",
        "E",
        "top_left",
    )
    s5e_panel_out = output_dir / "S5E_panel.pdf"
    panel_fig.savefig(s5e_panel_out, bbox_inches=None, facecolor="white")
    plt.close(panel_fig)

    results["S5_n_spots"] = float(len(data))
    results["S5_n_responder_spots"] = float(np.sum(data["response"] == "Responder"))
    results["S5_n_nonresponder_spots"] = float(np.sum(data["response"] == "NonResponder"))
    results["S5_file"] = inspect_tiff(s5_path)
    return results


def inspect_tiff(path: Path) -> dict[str, object]:
    with Image.open(path) as image:
        raw_dpi = image.info.get("dpi", (None, None))
        dpi = tuple(float(value) if value is not None else None for value in raw_dpi)
        try:
            display_path = str(path.relative_to(repository_root()))
        except ValueError:
            display_path = str(path)
        return {
            "sha256": sha256(path),
            "path": display_path,
            "size_pixels": list(image.size),
            "dpi": list(dpi),
            "size_inches": [image.size[0] / dpi[0], image.size[1] / dpi[1]],
            "mode": image.mode,
            "compression": image.info.get("compression"),
            "file_size_bytes": path.stat().st_size,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Override submission_figures")
    args = parser.parse_args()

    root = repository_root()
    data_dir = root / "processed_data"
    output_dir = args.output_dir.resolve() if args.output_dir else root / "submission_figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    configure_style()
    results: dict[str, object] = {}
    results.update(build_fig5(data_dir, output_dir))
    results.update(build_s5(data_dir, output_dir))
    results["declared_font_sizes_pt"] = FONT_SIZES
    results["environment"] = {
        "python": platform.python_version(),
        "matplotlib": mpl.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "seaborn": sns.__version__,
        "font_family": mpl.rcParams["font.family"],
    }
    manifest_path = output_dir / "fig5_s5_build_manifest.json"
    manifest_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
