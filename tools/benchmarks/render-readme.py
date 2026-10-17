#!/usr/bin/env python3
"""Render README tables and vector figures from benchmark-readme.R results.

Requires Python 3, NumPy, and Matplotlib. It never runs a benchmark or changes
measurements. Example:
  python tools/benchmarks/render-readme.py --results tools/benchmarks/results \
    --figures man/figures --pdfs tools/benchmarks/figures \
    --previews /tmp/rsvd-previews --update-readme README.md
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, MaxNLocator, NullFormatter
import numpy as np


COLORS = {"deterministic": "#263C53", "randomized": "#087F8C", "dense": "#BD5B29"}
STYLES = {
    "base_svd": ("base::svd", COLORS["deterministic"], "o", "-"),
    "rsvd": ("rsvd", COLORS["randomized"], "s", "-"),
    "prcomp": ("prcomp (dense)", COLORS["deterministic"], "o", "-"),
    "rpca_dense": ("rpca after densifying", COLORS["dense"], "^", "--"),
    "rpca": ("rpca, implicit sparse", COLORS["randomized"], "s", "-"),
}


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def short_number(value: float, _position: float = 0) -> str:
    if value >= 1000:
        return f"{value / 1000:g}k"
    if value >= 1:
        return f"{value:g}"
    return f"{value:.3g}"


def style_axes(
    ax: plt.Axes,
    sizes: list[int],
    ylabel: str,
    column_ratio: float,
    *,
    ratio: bool = False,
    xlabel: bool = True,
) -> None:
    ax.set_xscale("log", base=2)
    if not ratio:
        ax.set_yscale("log", base=10)
    ax.xaxis.set_major_locator(FixedLocator(sizes))
    ax.xaxis.set_major_formatter(FuncFormatter(short_number))
    ax.xaxis.set_minor_formatter(NullFormatter())
    if ratio:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}×"))
    else:
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=9))
        ax.yaxis.set_major_formatter(FuncFormatter(short_number))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", which="major", color="#DCE2E7", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#99A7B1")
    ax.spines[["left", "bottom"]].set_linewidth(0.75)
    ax.tick_params(axis="both", which="major", length=3.5, color="#75838E")
    ax.tick_params(axis="both", which="minor", length=2, color="#ABB5BD")
    if xlabel:
        ax.set_xlabel(f"Rows, m   (columns = {column_ratio:g} m)", labelpad=8)
    ax.set_ylabel(ylabel, labelpad=9)
    ax.set_xlim(min(sizes) / 1.08, max(sizes) * 1.20)


def ratio_limits(ax: plt.Axes, values: list[float]) -> None:
    upper = 1.20 * max(1.0, max(values))
    ax.set_ylim(0, upper)
    ticks = MaxNLocator(nbins=5, steps=(1, 2, 2.5, 5, 10)).tick_values(0, upper)
    ax.set_yticks([1.0] + [tick for tick in ticks if 1 < tick < upper])
    ax.axhline(1, color=COLORS["deterministic"], linestyle=":", linewidth=1.1,
               alpha=0.85, zorder=1)


def scaling_figure(
    rows: list[dict[str, str]],
    family: str,
    methods: list[str],
    title: str,
    subtitle: str,
    stem: str,
    args: argparse.Namespace,
) -> None:
    selected = [r for r in rows if r["family"] == family and r["scaling"] == "TRUE"]
    sizes = sorted({int(r["m"]) for r in selected})
    if len(sizes) < 3:
        raise ValueError(f"Need at least three measured sizes for {family}")
    repetitions = {int(r["repetitions"]) for r in selected}
    if len(repetitions) != 1:
        raise ValueError(f"Inconsistent repetition counts for {family}")
    reps = repetitions.pop()
    settings = [{int(r[key]) for r in selected} for key in ("k", "p", "q")]
    ratios = {int(r["n"]) / int(r["m"]) for r in selected}
    if any(len(values) != 1 for values in settings) or len(ratios) != 1:
        raise ValueError(f"Inconsistent rank, iteration, or aspect-ratio settings for {family}")
    k, p, q = (values.pop() for values in settings)
    column_ratio = ratios.pop()
    entries = [int(r["m"]) * int(r["n"]) / 1e6 for r in selected]
    nonzeros = [number(r, "input_nonzeros") / 1e6 for r in selected]
    densities = {number(r, "density") for r in selected}
    density_label = f"{100 * densities.pop():g}%" if len(densities) == 1 else "varying density"
    subtitle = subtitle.format(
        entries=f"{min(entries):g}–{max(entries):g}",
        nonzeros=f"{min(nonzeros):g}–{max(nonzeros):g}",
        density=density_label, k=k,
    )
    threads = args.metadata.get("OPENBLAS_NUM_THREADS", "unspecified")
    reference = methods[0]
    reference_label = {"base_svd": "base::svd", "prcomp": "prcomp"}[reference]
    measurements = {}
    for method in methods:
        data = sorted((r for r in selected if r["method"] == method), key=lambda r: int(r["m"]))
        if [int(r["m"]) for r in data] != sizes:
            raise ValueError(f"Incomplete size grid for {method}")
        values = {name: np.array([number(row, field) for row in data]) for name, field in (
            ("elapsed", "elapsed_median_s"),
            ("lower", "elapsed_q25_s"), ("upper", "elapsed_q75_s"),
            ("minimum", "elapsed_min_s"), ("maximum", "elapsed_max_s"),
            ("allocation_bytes", "r_allocation_total_bytes"),
        )}
        if any(not np.isfinite(value).all() or value.min() <= 0 for value in values.values()):
            raise ValueError("Plots and ratios require finite positive measured values")
        if (np.any(values["minimum"] > values["lower"]) or
            np.any(values["lower"] > values["elapsed"]) or
            np.any(values["elapsed"] > values["upper"]) or
            np.any(values["upper"] > values["maximum"])):
            raise ValueError("Inconsistent timing quantiles or extrema")
        measurements[method] = values

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.0))
    time_ax, speedup_ax = axes[0]
    memory_ax, reduction_ax = axes[1]
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.215, top=0.82,
                        wspace=0.31, hspace=0.43)
    fig.suptitle(title, x=0.095, y=0.985, ha="left", fontsize=17, fontweight="semibold")
    fig.text(0.095, 0.943, subtitle, ha="left", fontsize=11.2, color="#465866")
    handles = []
    time_limits = []
    memory_limits = []
    speedups = []
    reductions = []
    x = np.array(sizes)
    baseline = measurements[reference]
    for method in methods:
        values = measurements[method]
        elapsed = values["elapsed"]
        memory = values["allocation_bytes"] / 2**30
        label, color, marker, line = STYLES[method]
        time_ax.fill_between(x, values["lower"], values["upper"],
                             color=color, alpha=0.13, linewidth=0, zorder=1)
        time_ax.errorbar(x, elapsed,
                         yerr=(elapsed - values["minimum"], values["maximum"] - elapsed),
                         fmt="none", ecolor=color, elinewidth=0.9, capsize=3,
                         alpha=0.55, zorder=2)
        for ax, y in ((time_ax, elapsed), (memory_ax, memory)):
            ax.plot(x, y, color=color, marker=marker, linestyle=line, linewidth=2,
                    markersize=6, markeredgewidth=0.7, markeredgecolor="white", zorder=3)
        handles.append(Line2D([0], [0], color=color, marker=marker, linestyle=line,
                              linewidth=2, markersize=6, label=label))
        time_limits.extend([values["minimum"].min(), values["maximum"].max()])
        memory_limits.extend([memory.min(), memory.max()])

        if method != reference:
            # Ratios use the saved, unrounded medians or allocation totals.
            # They are not medians of paired ratios, and receive no IQR band.
            speedup = baseline["elapsed"] / elapsed
            reduction = baseline["allocation_bytes"] / values["allocation_bytes"]
            speedups.extend(speedup)
            reductions.extend(reduction)
            for ax, y in ((speedup_ax, speedup), (reduction_ax, reduction)):
                ax.plot(x, y, color=color, marker=marker, linestyle=line, linewidth=2,
                        markersize=6, markeredgewidth=0.7, markeredgecolor="white", zorder=3)
                ax.annotate(f"{y[-1]:.1f}×", (x[-1], y[-1]), xytext=(8, 5),
                            textcoords="offset points", ha="left", va="center",
                            fontsize=11, fontweight="semibold", color=color)

    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.084, 0.912),
               ncol=len(methods), frameon=False, columnspacing=2, handlelength=2.2,
               fontsize=11)
    style_axes(time_ax, sizes, "Seconds per call", column_ratio, xlabel=False)
    style_axes(speedup_ax, sizes, "Speedup (×)", column_ratio, ratio=True, xlabel=False)
    style_axes(memory_ax, sizes, "R allocations per call (GiB)", column_ratio)
    style_axes(reduction_ax, sizes, "Allocation reduction (×)", column_ratio, ratio=True)
    for ax, heading in (
        (time_ax, "A   Elapsed time"), (speedup_ax, "B   Speedup"),
        (memory_ax, "C   Cumulative allocation"), (reduction_ax, "D   Allocation reduction"),
    ):
        ax.set_title(heading, loc="left", fontsize=12, fontweight="semibold", pad=12)
    time_ax.set_ylim(min(time_limits) / 1.45, max(time_limits) * 1.45)
    memory_ax.set_ylim(min(memory_limits) / 1.45, max(memory_limits) * 1.45)
    ratio_limits(speedup_ax, speedups)
    ratio_limits(reduction_ax, reductions)
    randomized = [r for r in selected if r["method"] in ("rsvd", "rpca", "rpca_dense")]
    max_ratio = max(number(r, "error_over_optimal_max") for r in randomized)
    notes = (
        f"k = {k}, p = {p}, q = {q}  ·  {reps} calls per method  ·  {threads} BLAS threads  ·  reference: {reference_label}",
        "Time: medians, IQR bands, min–max whiskers. Speedup is the ratio of timing medians.",
        "Ratios = reference / method; 1× means equal cost. Allocations: one R profile each, not peak RSS.",
        f"Largest observed residual / optimal rank-{k} residual: {max_ratio:.6f}. Setup is excluded.",
    )
    for y, note in zip((0.126, 0.095, 0.064, 0.033), notes):
        fig.text(0.095, y, note, fontsize=9.5, color="#465866")
    args.figures.mkdir(parents=True, exist_ok=True)
    args.pdfs.mkdir(parents=True, exist_ok=True)
    args.previews.mkdir(parents=True, exist_ok=True)
    svg = args.figures / f"{stem}.svg"
    fig.savefig(svg, metadata={"Title": title})
    # Matplotlib leaves trailing spaces in multiline SVG path attributes.
    svg.write_text("\n".join(line.rstrip() for line in
                            svg.read_text(encoding="utf-8").splitlines()) + "\n",
                   encoding="utf-8")
    fig.savefig(args.pdfs / f"{stem}.pdf", metadata={"Title": title, "Creator": "rsvd benchmark renderer"})
    fig.savefig(args.previews / f"{stem}.png", dpi=220)
    plt.close(fig)


def value_pair(first: dict[str, str], second: dict[str, str], field: str, kind: str) -> str:
    values = [number(r, field) for r in (first, second)]
    if kind == "seconds":
        parts = [f"{v:.3f}" for v in values]
    elif kind == "gib":
        parts = [f"{v / 2**30:,.3f}" for v in values]
    elif kind == "percent":
        parts = [f"{100 * v:.3g}%" for v in values]
    else:
        raise ValueError(kind)
    return " / ".join(parts)


def benchmark_table(index: dict, entries: list[tuple[str, str, str, str]]) -> str:
    lines = [
        "| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for case, reference, method, label in entries:
        first, second = index[(case, reference)], index[(case, method)]
        m, n, k = (int(first[name]) for name in ("m", "n", "k"))
        lines.append(
            f"| {label}, {m:,} × {n:,} | {k} | "
            f"{value_pair(first, second, 'elapsed_median_s', 'seconds')} | "
            f"**{number(first, 'elapsed_median_s') / number(second, 'elapsed_median_s'):.1f}×** | "
            f"{value_pair(first, second, 'r_allocation_total_bytes', 'gib')} | "
            f"{value_pair(first, second, 'relative_error_median', 'percent')} |"
        )
    return "\n".join(lines)


def tables(rows: list[dict[str, str]]) -> dict[str, str]:
    index = {(r["case"], r["method"]): r for r in rows}
    return {
        "SVD": benchmark_table(index, [
            ("svd_tall_10000", "base_svd", "rsvd", "Dense"),
            ("svd_tall_20000", "base_svd", "rsvd", "Dense"),
        ]),
        "PCA": benchmark_table(index, [
            ("pca_dense", "prcomp", "rpca", "Dense PCA"),
            ("pca_sparse_scale_10000", "prcomp", "rpca", "Sparse-input PCA"),
            ("pca_sparse_scale_20000", "prcomp", "rpca", "Sparse-input PCA"),
        ]),
        "ID_CUR": benchmark_table(index, [
            ("id_dense", "id_deterministic", "rid", "Column ID, dense"),
            ("cur_dense", "cur_deterministic", "rcur", "CUR, dense"),
            ("id_sparse", "id_deterministic", "rid", "Column ID, sparse"),
            ("cur_sparse", "cur_deterministic", "rcur", "CUR, sparse"),
        ]),
    }


def replace_section(readme: str, name: str, body: str) -> str:
    start = f"<!-- BEGIN BENCHMARK {name} -->"
    end = f"<!-- END BENCHMARK {name} -->"
    if readme.count(start) != 1 or readme.count(end) != 1:
        raise ValueError(f"Expected exactly one marker pair for {name}")
    before, remaining = readme.split(start, 1)
    _old, after = remaining.split(end, 1)
    return before + start + "\n\n" + body + "\n\n" + end + after


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--figures", type=Path, default=Path("man/figures"))
    parser.add_argument("--pdfs", type=Path, default=Path("tools/benchmarks/figures"))
    parser.add_argument("--previews", type=Path, required=True)
    parser.add_argument("--update-readme", type=Path)
    args = parser.parse_args()
    with (args.results / "summary.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    with (args.results / "metadata.csv").open(newline="") as source:
        args.metadata = {row["key"]: row["value"] for row in csv.DictReader(source)}
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelsize": 11.5, "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
        "text.color": "#182C3B", "axes.labelcolor": "#263C53",
        "xtick.color": "#465866", "ytick.color": "#465866",
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    scaling_figure(rows, "svd", ["base_svd", "rsvd"],
                   "Dense SVD as matrix dimensions grow",
                   "{entries} million entries; geometric low-rank signal plus Gaussian noise",
                   "benchmark-svd-scaling", args)
    scaling_figure(rows, "pca", ["prcomp", "rpca_dense", "rpca"],
                   "PCA on sparse input: preserving sparse storage",
                   "{nonzeros} million nonzeros ({density}); same preprocessing and {k} score columns",
                   "benchmark-pca-scaling", args)
    blocks = tables(rows)
    legend = (
        "Paired values are **deterministic / randomized**. **Speedup** is the ratio of "
        "the deterministic and randomized median elapsed times, computed from unrounded "
        "measurements. R allocations are cumulative allocation volume in GiB, not peak RSS."
    )
    rendered = legend + "\n\n" + "\n\n".join(
        f"## {name}\n\n{body}" for name, body in blocks.items()) + "\n"
    (args.results / "tables.md").write_text(rendered, encoding="utf-8")
    (args.results / "figure-environment.json").write_text(json.dumps({
        "matplotlib": matplotlib.__version__, "numpy": np.__version__,
        "font": "DejaVu Sans", "svg_text": "editable text", "pdf_fonts": "embedded TrueType",
        "timing_band": "25th to 75th percentile of the recorded repetitions",
        "timing_whiskers": "minimum to maximum of all recorded repetitions",
        "repetition_counts": sorted({int(r["repetitions"]) for r in rows}),
        "allocation_unit": "GiB (2^30 bytes)",
        "blas_thread_limit": args.metadata.get("OPENBLAS_NUM_THREADS", "unspecified"),
        "figure_layout": "2 by 2: elapsed time, speedup, cumulative allocations, allocation reduction",
        "axis_scales": {"rows": "logarithmic base 2", "absolute_costs": "logarithmic base 10", "ratios": "linear"},
        "speedup": "deterministic median elapsed time / method median elapsed time; unrounded source values",
        "allocation_reduction": "deterministic R allocation bytes / method R allocation bytes; one profile each",
        "ratio_intervals": "none; timing quantile bands are not propagated into ratios",
        "ratio_reference": {"svd": "base::svd", "pca": "prcomp(as.matrix(A), ...)"},
    }, indent=2) + "\n", encoding="utf-8")
    if args.update_readme:
        readme = args.update_readme.read_text(encoding="utf-8")
        for name, block in blocks.items():
            readme = replace_section(readme, name, block)
        args.update_readme.write_text(readme, encoding="utf-8")
    print(f"Rendered {len(blocks)} tables and two SVG/PDF figures from {len(rows)} measured rows.")


if __name__ == "__main__":
    main()
