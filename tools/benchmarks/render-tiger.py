#!/usr/bin/env python3
"""Render the bundled tiger and measured low-rank reconstructions.

Run illustrate-tiger.R first. Numeric intermediates are kept outside the
repository; the small CSVs record the measurements and figure provenance.
Requires NumPy, Matplotlib, and Pillow (Matplotlib's PNG dependency).

  python tools/benchmarks/render-tiger.py \
    --results tools/benchmarks/tiger-results --images /tmp/rsvd-tiger \
    --figures man/figures --pdfs tools/benchmarks/figures
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


TITLE = "More rank, more detail"
CREDIT = "S. Taheri (2006), edited by Fir0002 · CC BY-SA 2.5 · Grayscale and low-rank adaptations"
SOURCE_URL = "https://commons.wikimedia.org/wiki/File:Sibirischer_tiger_de_edit02.jpg"
LICENSE_URL = "https://creativecommons.org/licenses/by-sa/2.5/"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--images", type=Path, required=True)
    parser.add_argument("--figures", type=Path, default=Path("man/figures"))
    parser.add_argument("--pdfs", type=Path, default=Path("tools/benchmarks/figures"))
    args = parser.parse_args()
    with (args.results / "summary.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    if [row["panel"] for row in rows] != ["original", "rank_10", "rank_50", "rank_150"]:
        raise ValueError("Expected original and measured ranks 10, 50, 150")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "text.color": "#263C53", "axes.titlecolor": "#263C53",
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    # Two rows keep the fine fur and whiskers readable at ordinary README width.
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.1))
    fig.subplots_adjust(left=0.045, right=0.955, top=0.875, bottom=0.12,
                        wspace=0.075, hspace=0.265)
    fig.suptitle(TITLE, x=0.045, y=0.982, ha="left", fontsize=17,
                 fontweight="semibold")
    fig.text(0.045, 0.943,
             "The same 1,600 × 1,200 image, approximated with fewer numbers.",
             ha="left", fontsize=11, color="#465866")
    for ax, row in zip(axes.flat, rows):
        m, n = int(row["m"]), int(row["n"])
        path = args.images / row["numeric_file"]
        numeric = np.fromfile(path, dtype="<f8")
        if numeric.size != m * n or not np.isfinite(numeric).all():
            raise ValueError(f"Invalid numeric image: {path}")
        matrix = numeric.reshape((m, n), order="F")
        fraction = float(row["factor_entry_fraction"])
        error = float(row["relative_frobenius_error"])
        if row["panel"] == "original":
            title = "Original"
            annotation = "1,920,000 image entries · reference"
        else:
            k = int(row["k"])
            expected = k * (m + n + 1) / (m * n)
            if not np.isclose(fraction, expected, rtol=1e-12):
                raise ValueError("Stored factor entry fraction is inconsistent")
            title = f"Rank {k}"
            annotation = f"{100 * fraction:.2f}% of entries · {100 * error:.1f}% error"
        ax.imshow(matrix.T, cmap="gray", vmin=0, vmax=1, origin="lower",
                  interpolation="lanczos", aspect="equal", rasterized=True)
        ax.set_axis_off()
        ax.set_title(title, loc="left", fontsize=12, fontweight="semibold", pad=9)
        ax.text(0.0, -0.057, annotation, transform=ax.transAxes, va="top",
                ha="left", fontsize=10.5,
                color="#465866" if row["panel"] == "original" else "#087F8C")

    fig.text(0.045, 0.066,
             "rsvd: p = 10, q = 2, seed = 42 for each rank. Error: relative Frobenius norm before display clipping.",
             fontsize=8.7, color="#465866")
    fig.text(0.045, 0.043,
             "Factor storage counts k(m + n + 1) numeric entries; these percentages are not JPEG or file-compression ratios.",
             fontsize=8.7, color="#465866")
    fig.text(0.045, 0.020, CREDIT, fontsize=8.0, color="#687782")
    args.figures.mkdir(parents=True, exist_ok=True)
    args.pdfs.mkdir(parents=True, exist_ok=True)
    png = args.figures / "tiger-rank-comparison.png"
    pdf = args.pdfs / "tiger-rank-comparison.pdf"
    description = (
        "Actual bundled tiger data and independent rsvd rank-10, rank-50, and "
        "rank-150 reconstructions. " + CREDIT + ". Source: " + SOURCE_URL +
        ". License: " + LICENSE_URL + "."
    )
    # A README-sized preview avoids bloating the installed help/figures tree;
    # the separate PDF retains the higher-resolution image export.
    fig.savefig(png, dpi=120, metadata={"Title": TITLE, "Description": description},
                pil_kwargs={"optimize": True})
    fig.savefig(pdf, dpi=180, metadata={
        "Title": TITLE, "Author": "rsvd contributors; image by S. Taheri, edited by Fir0002",
        "Subject": description, "Creator": "rsvd render-tiger.py",
    })
    plt.close(fig)
    with Image.open(png) as rendered:
        png_size = list(rendered.size)
    environment = {
        "python": platform.python_version(), "numpy": np.__version__,
        "matplotlib": matplotlib.__version__, "png_pixels": png_size,
        "display_limits": [0, 1], "display_origin": "lower",
        "display_orientation": "transpose of R matrix, equivalent to image(tiger)",
        "raster_interpolation": "lanczos; same for all panels",
        "png_dpi": 120, "pdf_dpi": 180, "pdf_fonts": "TrueType",
        "png_sha256": hashlib.sha256(png.read_bytes()).hexdigest(),
        "source_url": SOURCE_URL, "license_url": LICENSE_URL,
        "image_credit": CREDIT,
    }
    (args.results / "figure-environment.json").write_text(
        json.dumps(environment, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
