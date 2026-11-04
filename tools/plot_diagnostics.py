#!/usr/bin/env python3
"""Plot front history and conservation diagnostics from output_index.csv."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_history(path):
    columns = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        for name in reader.fieldnames:
            columns[name] = []
        for row in reader:
            for name in reader.fieldnames:
                columns[name].append(row[name])
    numeric = {
        "time_us": np.asarray(columns["time_us"], dtype=float),
        "step": np.asarray(columns["step"], dtype=int),
        "total_cells": np.asarray(columns["total_cells"], dtype=float),
        "front_x_cm": np.asarray(columns["front_x_m"], dtype=float) * 100.0,
        "reaction_front_x_cm": (
            np.asarray(columns["reaction_front_x_m"], dtype=float) * 100.0),
        "separation_mm": np.asarray(
            columns["front_separation_mm"], dtype=float),
        "mass": np.asarray(columns["mass"], dtype=float),
        "energy": np.asarray(columns["energy"], dtype=float),
    }
    return numeric


def gradient_speed(time_us, position_cm):
    t = time_us * 1.0e-6
    x = position_cm * 1.0e-2
    return np.gradient(x, t)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_index", type=Path)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()

    path = args.output_index.resolve()
    history = load_history(path)
    outdir = args.outdir.resolve() if args.outdir else path.parent / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    })

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                             constrained_layout=True)
    axes[0].plot(
        history["time_us"], history["front_x_cm"],
        "o-", label="pressure front")
    axes[0].plot(
        history["time_us"], history["reaction_front_x_cm"],
        "s-", label="reaction front")
    axes[0].set_ylabel("front position (cm)")
    axes[0].grid(alpha=0.22)
    axes[0].legend()

    shock_speed = gradient_speed(
        history["time_us"], history["front_x_cm"])
    reaction_speed = gradient_speed(
        history["time_us"], history["reaction_front_x_cm"])
    axes[1].plot(
        history["time_us"], shock_speed, "o-",
        label="pressure-front speed")
    axes[1].plot(
        history["time_us"], reaction_speed, "s-",
        label="reaction-front speed")
    axes[1].set_xlabel(r"time ($\mu$s)")
    axes[1].set_ylabel("speed (m/s)")
    axes[1].grid(alpha=0.22)
    axes[1].legend()
    fig.suptitle("Detonation front history")
    fig.savefig(outdir / "front_history.png", dpi=200)
    plt.close(fig)

    mass_drift = history["mass"] / history["mass"][0] - 1.0
    energy_drift = history["energy"] / history["energy"][0] - 1.0
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                             constrained_layout=True)
    axes[0].plot(
        history["time_us"], history["separation_mm"], "o-", color="tab:red")
    axes[0].set_ylabel("front separation (mm)")
    axes[0].grid(alpha=0.22)
    axes[1].plot(history["time_us"], mass_drift, "o-", label="mass")
    axes[1].plot(history["time_us"], energy_drift, "s-", label="energy")
    axes[1].set_ylabel("relative drift")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[1].grid(alpha=0.22)
    axes[1].legend()
    axes[2].plot(
        history["time_us"], history["total_cells"] / 1.0e6,
        "o-", color="tab:green")
    axes[2].set_xlabel(r"time ($\mu$s)")
    axes[2].set_ylabel("AMR cells (million)")
    axes[2].grid(alpha=0.22)
    fig.suptitle("Coupling, conservation and AMR history")
    fig.savefig(outdir / "conservation_and_separation.png", dpi=200)
    plt.close(fig)

    print("time_us  shock_m_s  reaction_m_s  separation_mm")
    for values in zip(
            history["time_us"], shock_speed, reaction_speed,
            history["separation_mm"]):
        print("%8.3f %11.1f %13.1f %15.3f" % values)
    print("mass_drift", mass_drift[-1])
    print("energy_drift", energy_drift[-1])
    print("Wrote", outdir / "front_history.png")
    print("Wrote", outdir / "conservation_and_separation.png")


if __name__ == "__main__":
    main()
