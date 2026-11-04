#!/usr/bin/env python3
"""Plot AMR box layout and refinement level from an AMReX plotfile."""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle


BOX_PATTERN = re.compile(
    r"\(\((\d+),(\d+)\)\s+\((\d+),(\d+)\)\s+\(0,0\)\)")


def read_domain(plotdir: Path):
    lines = (plotdir / "Header").read_text(encoding="utf-8").splitlines()
    nvars = int(lines[1].strip())
    location = 2 + nvars
    time_us = float(lines[location + 1]) * 1e6
    problo = [float(value) for value in lines[location + 3].split()]
    probhi = [float(value) for value in lines[location + 4].split()]
    levels = int(lines[location + 2]) + 1
    ratios = [int(value) for value in lines[location + 5].split()]
    spacing = [
        [float(value) for value in lines[location + 8 + lev].split()]
        for lev in range(levels)
    ]
    return time_us, problo, probhi, levels, ratios, spacing


def read_boxes(header: Path):
    text = header.read_text(encoding="utf-8")
    return [
        tuple(int(value) for value in match)
        for match in BOX_PATTERN.findall(text)
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--xlim", type=float, default=0.06)
    args = parser.parse_args()
    plotdir = args.plotfile.resolve()
    outdir = plotdir.parent / "figures"
    outdir.mkdir(exist_ok=True)
    time_us, problo, probhi, levels, ratios, spacing = read_domain(plotdir)
    level_boxes = [
        read_boxes(plotdir / ("Level_%d" % lev) / "Cell_H")
        for lev in range(levels)
    ]
    level_cells = [
        sum((i1 - i0 + 1) * (j1 - j0 + 1) for i0, j0, i1, j1 in boxes)
        for boxes in level_boxes
    ]

    colors = ["#4c78a8", "#59a14f", "#f28e2b", "#e15759", "#b07aa1"]
    labels = [
        "Level %d: %d cells" % (lev, level_cells[lev])
        for lev in range(levels)
    ]
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    })

    fig, axes = plt.subplots(
        2, 1, figsize=(13, 7.6), constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 1.15]})
    for ax, xmax in zip(axes, [probhi[0], args.xlim]):
        for lev, boxes in enumerate(level_boxes):
            dx = spacing[lev][0]
            dy = spacing[lev][1]
            for i0, j0, i1, j1 in boxes:
                x0 = problo[0] + i0 * dx
                x1 = problo[0] + (i1 + 1) * dx
                y0 = problo[1] + j0 * dy
                y1 = problo[1] + (j1 + 1) * dy
                ax.add_patch(Rectangle(
                    (x0, y0 * 1e3), x1 - x0, (y1 - y0) * 1e3,
                    facecolor=colors[lev], edgecolor="white",
                    linewidth=0.45, alpha=0.75, zorder=lev + 1))
        ax.set_xlim(problo[0], xmax)
        ax.set_ylim(problo[1] * 1e3, probhi[1] * 1e3)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (mm)")
        ax.grid(color="0.85", linewidth=0.35, alpha=0.45)
    axes[0].set_title("AMR coverage over the full domain")
    axes[1].set_title("Wave-front refinement detail")
    axes[0].legend(
        handles=[
            Patch(facecolor=colors[lev], edgecolor="white", alpha=0.75,
                  label=labels[lev])
            for lev in range(levels)
        ],
        loc="upper right", ncol=2, framealpha=0.92)
    fig.suptitle(
        "AMReX grid layout at t = %.3f $\\mu$s" % time_us)
    output = outdir / "amr_levels.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print("Wrote", output)
    print("level cells:", level_cells)


if __name__ == "__main__":
    main()
