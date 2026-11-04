#!/usr/bin/env python3
"""Create final flow-structure and centerline figures from saved field data."""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


BOX_PATTERN = re.compile(
    r"\(\((\d+),(\d+)\)\s+\((\d+),(\d+)\)\s+\(0,0\)\)")


def read_header(plotdir: Path):
    lines = (plotdir / "Header").read_text(encoding="utf-8").splitlines()
    nvars = int(lines[1].strip())
    loc = nvars + 2
    time_us = float(lines[loc + 1]) * 1e6
    problo = np.fromstring(lines[loc + 3], sep=" ")
    probhi = np.fromstring(lines[loc + 4], sep=" ")
    nlevels = int(lines[loc + 2]) + 1
    spacing = [
        np.fromstring(lines[loc + 8 + lev], sep=" ")
        for lev in range(nlevels)
    ]
    return time_us, problo, probhi, nlevels, spacing


def read_boxes(path: Path):
    return [
        tuple(int(v) for v in match)
        for match in BOX_PATTERN.findall(path.read_text(encoding="utf-8"))
    ]


def level_map(plotdir, nlevels, problo, probhi, spacing, x, y):
    result = np.zeros((y.size, x.size), dtype=np.int16)
    for lev in range(nlevels):
        dx = spacing[lev][0]
        dy = spacing[lev][1]
        for i0, j0, i1, j1 in read_boxes(
                plotdir / ("Level_%d" % lev) / "Cell_H"):
            x0 = problo[0] + i0 * dx
            x1 = problo[0] + (i1 + 1) * dx
            y0 = problo[1] + j0 * dy
            y1 = problo[1] + (j1 + 1) * dy
            ix0 = np.searchsorted(x, x0, side="left")
            ix1 = np.searchsorted(x, x1, side="right")
            iy0 = np.searchsorted(y, y0, side="left")
            iy1 = np.searchsorted(y, y1, side="right")
            result[iy0:iy1, ix0:ix1] = lev
    return result


def panel(ax, x, y, value, title, label, cmap, vmin=None, vmax=None):
    image = ax.imshow(
        value, origin="lower",
        extent=(x[0], x[-1], y[0] * 1e3, y[-1] * 1e3),
        aspect="auto", interpolation="nearest", cmap=cmap,
        vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (mm)")
    cb = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label(label)
    cb.ax.tick_params(labelsize=8)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument("--data", type=Path)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    plotdir = args.plotfile.resolve()
    outdir = (
        args.outdir.resolve() if args.outdir is not None
        else plotdir.parent / "figures")
    outdir.mkdir(parents=True, exist_ok=True)
    time_us, problo, probhi, nlevels, spacing = read_header(plotdir)
    data_path = args.data
    if data_path is None:
        data_path = outdir / ("fields_%.3fus.npz" % time_us)
    data = np.load(data_path)
    x = data["x"]
    y = data["y"]
    rho = data["rho"]
    p = data["p_atm"]
    temp = data["T"]
    y_oh = data["Y_OH"]
    y_o2 = data["Y_O2"]

    drho_dx = np.gradient(rho, x, axis=1)
    drho_dy = np.gradient(rho, y, axis=0)
    grad_rho = np.log10(
        1.0 + np.hypot(drho_dx, drho_dy) /
        np.maximum(rho, 1.0e-12))
    levels = level_map(plotdir, nlevels, problo, probhi, spacing, x, y)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    })
    fig = plt.figure(figsize=(16, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 3)
    ax_t = fig.add_subplot(grid[0, 0])
    ax_l = fig.add_subplot(grid[0, 1])
    ax_p = fig.add_subplot(grid[0, 2])
    ax_g = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1:])

    panel(
        ax_t, x, y, temp, "Temperature", "T (K)", "inferno",
        np.percentile(temp, 0.2), np.percentile(temp, 99.8))
    level_colors = ListedColormap(
        ["#4c78a8", "#59a14f", "#f28e2b", "#e15759"])
    image = ax_l.imshow(
        levels, origin="lower",
        extent=(x[0], x[-1], y[0] * 1e3, y[-1] * 1e3),
        aspect="auto", interpolation="nearest",
        cmap=level_colors, norm=BoundaryNorm(
            np.arange(-0.5, nlevels + 0.5), nlevels))
    ax_l.set_title("AMR refinement level")
    ax_l.set_xlabel("x (m)")
    ax_l.set_ylabel("y (mm)")
    cb = plt.colorbar(
        image, ax=ax_l, fraction=0.046, pad=0.025,
        ticks=np.arange(nlevels))
    cb.set_label("level")
    panel(
        ax_p, x, y, p, "Pressure", "p (atm)", "magma",
        vmax=np.percentile(p, 99.8))
    panel(
        ax_g, x, y, grad_rho, "Normalized density gradient",
        r"$\log_{10}(1+|\nabla\rho|/\rho)$", "gray",
        np.percentile(grad_rho, 1), np.percentile(grad_rho, 99.7))

    xz = (x >= 0.05) & (x <= 0.13)
    xc = x[xz]
    panel(
        ax_c, xc, y, grad_rho[:, xz],
        "Cellular reaction-front structure",
        r"$\log_{10}(1+|\nabla\rho|/\rho)$", "gray",
        np.percentile(grad_rho, 1), np.percentile(grad_rho, 99.7))
    cs = ax_c.contour(
        xc, y * 1e3, y_oh[:, xz],
        levels=[0.01, 0.05, 0.10], colors=["#e15759", "#ffbf00", "#ffffff"],
        linewidths=[0.9, 1.1, 1.4])
    ax_c.clabel(cs, inline=True, fontsize=7, fmt="%.2f")
    fig.suptitle(
        "2D detonation flow structure at t = %.3f $\\mu$s" % time_us)
    time_tag = f"{time_us:.0f}us"
    flow_output = outdir / ("flow_structure_" + time_tag + ".png")
    fig.savefig(flow_output, dpi=220)
    plt.close(fig)

    mid = y.size // 2
    p_line = p[mid]
    T_line = temp[mid]
    o2_line = y_o2[mid]
    reaction_x = x[np.argmax(np.abs(np.gradient(o2_line, x)))]
    pressure_x = x[np.argmax(np.abs(np.gradient(p_line, x)))]

    fig, axes = plt.subplots(
        3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    axes[0].plot(x, p_line, color="tab:blue")
    axes[0].set_ylabel("p (atm)")
    axes[0].set_title("P-x, T-x and O2 mass fraction-x at t = %.3f $\\mu$s"
                      % time_us)
    axes[1].plot(x, T_line, color="tab:red")
    axes[1].set_ylabel("T (K)")
    axes[2].plot(x, o2_line, color="tab:green")
    axes[2].set_ylabel(r"$Y_{O_2}$")
    axes[2].set_xlabel("x (m)")
    for ax in axes:
        ax.axvline(reaction_x, color="0.35", linestyle="--", linewidth=0.9)
        ax.axvline(pressure_x, color="0.6", linestyle=":", linewidth=0.9)
        ax.grid(alpha=0.22)
        ax.set_xlim(problo[0], probhi[0])
    axes[0].text(
        reaction_x, axes[0].get_ylim()[1] * 0.93,
        " reaction front", va="top", fontsize=8, color="0.25")
    axes[0].text(
        pressure_x, axes[0].get_ylim()[1] * 0.84,
        " pressure front", va="top", fontsize=8, color="0.45")
    profile_output = outdir / ("profiles_" + time_tag + ".png")
    fig.savefig(profile_output, dpi=220)
    plt.close(fig)

    print("Wrote", flow_output)
    print("Wrote", profile_output)
    print("reaction_x", reaction_x, "pressure_x", pressure_x)


if __name__ == "__main__":
    main()
