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
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter
import numpy as np


COLORS = {"deterministic": "#263C53", "randomized": "#087F8C", "dense": "#BD5B29"}
STYLES = {
    "base_svd": ("Base SVD", COLORS["deterministic"], "o", "-"),
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


def style_axes(ax: plt.Axes, sizes: list[int], ylabel: str, column_ratio: float) -> None:
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=10)
    ax.xaxis.set_major_locator(FixedLocator(sizes))
    ax.xaxis.set_major_formatter(FuncFormatter(short_number))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.yaxis.set_major_formatter(FuncFormatter(short_number))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", which="major", color="#DCE2E7", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#99A7B1")
    ax.spines[["left", "bottom"]].set_linewidth(0.75)
    ax.tick_params(axis="both", which="major", length=3.5, color="#75838E")
    ax.tick_params(axis="both", which="minor", length=2, color="#ABB5BD")
    ax.set_xlabel(f"Rows, m   (columns = {column_ratio:g} m)", labelpad=8)
    ax.set_ylabel(ylabel, labelpad=9)
    ax.set_xlim(min(sizes) / 1.13, max(sizes) * 1.13)


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
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9))
    fig.subplots_adjust(left=0.087, right=0.981, bottom=0.225, top=0.72, wspace=0.29)
    fig.suptitle(title, x=0.087, y=0.985, ha="left", fontsize=15, fontweight="semibold")
    fig.text(0.087, 0.900, subtitle, ha="left", fontsize=10.5, color="#465866")
    handles = []
    time_limits = []
    memory_limits = []
    for method in methods:
        data = sorted((r for r in selected if r["method"] == method), key=lambda r: int(r["m"]))
        if [int(r["m"]) for r in data] != sizes:
            raise ValueError(f"Incomplete size grid for {method}")
        x = np.array(sizes)
        elapsed = np.array([number(r, "elapsed_median_s") for r in data])
        lower = np.array([number(r, "elapsed_q25_s") for r in data])
        upper = np.array([number(r, "elapsed_q75_s") for r in data])
        memory = np.array([number(r, "r_allocation_total_bytes") / 2**30 for r in data])
        if min(elapsed.min(), lower.min(), memory.min()) <= 0:
            raise ValueError("Logarithmic plots require positive measured values")
        label, color, marker, line = STYLES[method]
        axes[0].fill_between(x, lower, upper, color=color, alpha=0.13, linewidth=0)
        for ax, y in zip(axes, (elapsed, memory)):
            ax.plot(x, y, color=color, marker=marker, linestyle=line, linewidth=2,
                    markersize=5.5, markeredgewidth=0.7, markeredgecolor="white")
        handles.append(Line2D([0], [0], color=color, marker=marker, linestyle=line,
                              linewidth=2, markersize=5.5, label=label))
        time_limits.extend([lower.min(), upper.max()])
        memory_limits.extend([memory.min(), memory.max()])
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.077, 0.855),
               ncol=len(methods), frameon=False, columnspacing=2, handlelength=2.2,
               fontsize=10)
    style_axes(axes[0], sizes, "Seconds per call", column_ratio)
    style_axes(axes[1], sizes, "R allocations per call (GiB)", column_ratio)
    axes[0].set_title("A   Elapsed time", loc="left", fontsize=11, fontweight="semibold", pad=11)
    axes[1].set_title("B   Cumulative allocation", loc="left", fontsize=11, fontweight="semibold", pad=11)
    axes[0].set_ylim(min(time_limits) / 1.5, max(time_limits) * 1.5)
    axes[1].set_ylim(min(memory_limits) / 1.5, max(memory_limits) * 1.5)
    randomized = [r for r in selected if r["method"] in ("rsvd", "rpca", "rpca_dense")]
    max_ratio = max(number(r, "error_over_optimal_max") for r in randomized)
    fig.text(0.087, 0.082,
             f"k = {k}, p = {p}, q = {q}  ·  {reps}-run timing medians, IQR shaded  ·  allocations: one profile  ·  BLAS threads: {threads}",
             fontsize=9, color="#465866")
    fig.text(0.087, 0.037,
             f"Largest observed error / optimal rank-{k} error: {max_ratio:.4f}.  "
             "R allocations exclude setup and do not measure peak RSS.",
             fontsize=9, color="#465866")
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
        "| Input | k | Time (s), deterministic / randomized | R allocations (GiB), deterministic / randomized | Relative error, deterministic / randomized |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for case, reference, method, label in entries:
        first, second = index[(case, reference)], index[(case, method)]
        m, n, k = (int(first[name]) for name in ("m", "n", "k"))
        lines.append(
            f"| {label}, {m:,} × {n:,} | {k} | "
            f"{value_pair(first, second, 'elapsed_median_s', 'seconds')} | "
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
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
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
    rendered = "\n\n".join(f"## {name}\n\n{body}" for name, body in blocks.items()) + "\n"
    (args.results / "tables.md").write_text(rendered, encoding="utf-8")
    (args.results / "figure-environment.json").write_text(json.dumps({
        "matplotlib": matplotlib.__version__, "numpy": np.__version__,
        "font": "DejaVu Sans", "svg_text": "editable text", "pdf_fonts": "embedded TrueType",
        "timing_band": "25th to 75th percentile of the recorded repetitions",
        "repetition_counts": sorted({int(r["repetitions"]) for r in rows}),
        "allocation_unit": "GiB (2^30 bytes)",
        "blas_thread_limit": args.metadata.get("OPENBLAS_NUM_THREADS", "unspecified"),
    }, indent=2) + "\n", encoding="utf-8")
    if args.update_readme:
        readme = args.update_readme.read_text(encoding="utf-8")
        for name, block in blocks.items():
            readme = replace_section(readme, name, block)
        args.update_readme.write_text(readme, encoding="utf-8")
    print(f"Rendered {len(blocks)} tables and two SVG/PDF figures from {len(rows)} measured rows.")


if __name__ == "__main__":
    main()
