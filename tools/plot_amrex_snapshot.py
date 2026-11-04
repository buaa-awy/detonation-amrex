#!/usr/bin/env python3
"""Build 2D field figures from an AMReX plotfile using fextract."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RU = 8.31446261815324
P0 = 101325.0


def wsl_path(path: Path) -> str:
    text = str(path.resolve()).replace("\\", "/")
    if len(text) > 2 and text[1:3] == ":/":
        return "/mnt/" + text[0].lower() + text[2:]
    return text


def load_mechanism(root: Path):
    import sys

    sys.path.insert(0, str(root / "src"))
    from det.mech import load_mechanism

    return load_mechanism(str(root / "mechanism" / "mech_9sp21r.yaml")).nb


def temperature_from_energy(nb, energy, Y):
    lo = np.full(energy.shape, 250.0)
    hi = np.full(energy.shape, 6500.0)
    for _ in range(70):
        T = 0.5 * (lo + hi)
        a = np.where(
            (T[:, None] <= nb.Tmid[None, :])[:, :, None],
            nb.A_low[None, :, :],
            nb.A_high[None, :, :],
        )
        h_rt = (
            a[:, :, 0]
            + T[:, None]
            * (
                a[:, :, 1] / 2.0
                + T[:, None]
                * (
                    a[:, :, 2] / 3.0
                    + T[:, None]
                    * (a[:, :, 3] / 4.0 + T[:, None] * a[:, :, 4] / 5.0)
                )
            )
            + a[:, :, 5] / T[:, None]
        )
        e = np.sum(
            Y * (RU * T[:, None] * (h_rt - 1.0) / nb.W[None, :]),
            axis=1,
        )
        above = e > energy
        hi = np.where(above, T, hi)
        lo = np.where(above, lo, T)
    return 0.5 * (lo + hi)


def read_header(plotdir: Path):
    lines = (plotdir / "Header").read_text(encoding="utf-8").splitlines()
    lines = [line.strip() for line in lines]
    nvars = int(lines[1])
    names = lines[2:2 + nvars]
    location = 2 + nvars
    time = float(lines[location + 1])
    problo = np.fromstring(lines[location + 3], sep=" ")
    probhi = np.fromstring(lines[location + 4], sep=" ")
    base_dx = np.fromstring(lines[location + 8], sep=" ")
    return names, time, problo, probhi, base_dx


def resolve_fextract(user_value: Optional[str]) -> str:
    if user_value:
        return str(user_value)
    if os.environ.get("FEXTRACT"):
        return os.environ["FEXTRACT"]
    candidates = []
    if os.environ.get("AMREX_ROOT"):
        candidates.append(
            Path(os.environ["AMREX_ROOT"]) /
            "Tools/Plotfile/fextract.gnu.MPI.ex")
    candidates.extend([
        Path.home() / "amrex/Tools/Plotfile/fextract.gnu.MPI.ex",
        Path.home() / "amrex-install-2d-cuda/bin/fextract.gnu.MPI.ex",
        Path.home() / "amrex-install-2d-cpu/bin/fextract.gnu.MPI.ex",
    ])
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise FileNotFoundError(
        "AMReX fextract was not found. Set --fextract or FEXTRACT.")


def extract_slice(fextract: str, plotdir: Path, y: float, outfile: Path):
    cmd = [
        "wsl.exe", "-d", "Ubuntu-20.04", "--exec", fextract,
        "-s", wsl_path(outfile), "-d", "0", "-y", f"{y:.16g}",
        wsl_path(plotdir),
    ]
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        raise RuntimeError(proc.stdout + "\n" + proc.stderr)
    return np.loadtxt(outfile, encoding="utf-8")


def fields_from_plot(nb, data, names):
    col = {name: i + 1 for i, name in enumerate(names)}
    species = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "N2"]
    if all(name in col for name in (
            "p", "T", "U_x", "U_y", "rho")) and all(
                "Y_" + name in col for name in species):
        rho = data[:, col["rho"]]
        u = data[:, col["U_x"]]
        v = data[:, col["U_y"]]
        T = data[:, col["T"]]
        Y = np.stack(
            [data[:, col["Y_" + name]] for name in species], axis=1)
        progress = (
            data[:, col["progress"]] if "progress" in col
            else np.clip(1.0 - Y[:, 0] / 0.1119067440, 0.0, 1.0))
        return {
            "rho": rho,
            "u": u,
            "v": v,
            "T": T,
            "p_atm": data[:, col["p"]] / P0,
            "Y_H2": Y[:, 0],
            "Y_O2": Y[:, 3],
            "Y_OH": Y[:, 4],
            "Y_H2O": Y[:, 5],
            "progress": progress,
        }

    # Read legacy plotfiles that stored conserved variables only.
    rho = data[:, col["density"]]
    u = data[:, col["momentum_x"]] / rho
    v = data[:, col["momentum_y"]] / rho
    kinetic = 0.5 * rho * (u * u + v * v)
    internal = (data[:, col["total_energy"]] - kinetic) / rho
    species = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "N2"]
    Y = np.stack(
        [data[:, col["rho_" + s]] / rho for s in species], axis=1)
    Y = np.maximum(Y, 0.0)
    Y /= np.maximum(Y.sum(axis=1, keepdims=True), 1e-30)
    T = temperature_from_energy(nb, internal, Y)
    Rmix = np.sum(Y * RU / nb.W[None, :], axis=1)
    p = rho * Rmix * T
    return {
        "rho": rho,
        "u": u,
        "v": v,
        "T": T,
        "p_atm": p / P0,
        "Y_H2": Y[:, 0],
        "Y_O2": Y[:, 3],
        "Y_OH": Y[:, 4],
        "Y_H2O": Y[:, 5],
        "progress": np.clip(1.0 - Y[:, 0] / 0.1119067440, 0.0, 1.0),
    }


def field_panel(ax, x, y, values, title, label, cmap, vmin=None, vmax=None):
    image = ax.imshow(
        values,
        origin="lower",
        extent=(x[0], x[-1], y[0], y[-1]),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (mm)")
    ax.set_ylim(y[0], y[-1])
    ax.tick_params(labelsize=9)
    cb = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label)
    cb.ax.tick_params(labelsize=8)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("plotfile", type=Path)
    parser.add_argument(
        "--fextract",
        default=None,
    )
    parser.add_argument("--ny", type=int, default=0)
    parser.add_argument("--nx", type=int, default=3001)
    parser.add_argument("--zoom-center", type=float, default=None)
    parser.add_argument("--zoom-width", type=float, default=0.06)
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    plotdir = args.plotfile.resolve()
    outdir = plotdir.parent / "figures"
    outdir.mkdir(exist_ok=True)
    names, time, problo, probhi, base_dx = read_header(plotdir)
    if args.ny <= 0:
        args.ny = max(
            1, int(round((probhi[1] - problo[1]) / base_dx[1])))
    fextract = resolve_fextract(args.fextract)
    nb = load_mechanism(root)

    cache = outdir / ("fields_%.3fus.npz" % (time * 1e6))
    if cache.exists():
        with np.load(cache) as loaded:
            arrays = {key: loaded[key] for key in loaded.files}
    else:
        xout = np.linspace(problo[0], probhi[0], args.nx)
        yout = problo[1] + (np.arange(args.ny) + 0.5) * (
            (probhi[1] - problo[1]) / args.ny)
        values = {key: [] for key in (
            "rho", "u", "v", "T", "p_atm", "Y_H2", "Y_O2",
            "Y_OH", "Y_H2O", "progress")}

        def extract_one(j, y, tmp):
            slice_file = Path(tmp) / ("slice_%04d.dat" % j)
            data = extract_slice(fextract, plotdir, float(y), slice_file)
            if data.ndim == 1:
                data = data[None, :]
            x = data[:, 0]
            order = np.argsort(x)
            x = x[order]
            data = data[order]
            fields = fields_from_plot(nb, data, names)
            return j, {
                key: np.interp(xout, x, fields[key]) for key in values
            }

        with tempfile.TemporaryDirectory(dir=str(outdir)) as tmp:
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                futures = [
                    pool.submit(extract_one, j, y, tmp)
                    for j, y in enumerate(yout)
                ]
                completed = 0
                for future in as_completed(futures):
                    j, result = future.result()
                    for key, value in result.items():
                        values[key].append((j, value))
                    completed += 1
                    print(
                        "extracted %d/%d" % (completed, args.ny),
                        flush=True)

        arrays = {
            key: np.asarray([value for _, value in sorted(vals)])
            for key, vals in values.items()
        }
        arrays.update({
            "x": xout,
            "y": yout,
            "time_us": np.array(time * 1e6),
        })
        np.savez_compressed(cache, **arrays)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    })

    fig, axes = plt.subplots(2, 2, figsize=(12, 6.2), constrained_layout=True)
    field_panel(
        axes[0, 0], arrays["x"], arrays["y"] * 1e3, arrays["rho"],
        "Density", r"$\rho$ (kg/m$^3$)", "viridis")
    field_panel(
        axes[0, 1], arrays["x"], arrays["y"] * 1e3, arrays["p_atm"],
        "Pressure", "p (atm)", "magma",
        vmax=np.percentile(arrays["p_atm"], 99.7))
    field_panel(
        axes[1, 0], arrays["x"], arrays["y"] * 1e3, arrays["T"],
        "Temperature", "T (K)", "inferno",
        vmin=np.percentile(arrays["T"], 0.3),
        vmax=np.percentile(arrays["T"], 99.7))
    field_panel(
        axes[1, 1], arrays["x"], arrays["y"] * 1e3, arrays["progress"],
        "H2 consumption", "reaction progress", "cividis", 0.0, 1.0)
    fig.suptitle(r"2D detonation fields at t = %.3f $\mu$s" % (time * 1e6))
    fig.savefig(outdir / "fields_density_pressure_temperature.png", dpi=220)
    plt.close(fig)

    mid = args.ny // 2
    if args.zoom_center is None:
        zoom_center = float(arrays["x"][np.argmax(arrays["p_atm"][mid])])
    else:
        zoom_center = args.zoom_center
    zoom_lo = max(problo[0], zoom_center - 0.5 * args.zoom_width)
    zoom_hi = min(probhi[0], zoom_center + 0.5 * args.zoom_width)
    zoom = (arrays["x"] >= zoom_lo) & (arrays["x"] <= zoom_hi)
    xz = arrays["x"][zoom]
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.2), constrained_layout=True)
    field_panel(
        axes[0, 0], xz, arrays["y"] * 1e3, arrays["rho"][:, zoom],
        "Density", r"$\rho$ (kg/m$^3$)", "viridis")
    field_panel(
        axes[0, 1], xz, arrays["y"] * 1e3, arrays["p_atm"][:, zoom],
        "Pressure", "p (atm)", "magma",
        vmax=np.percentile(arrays["p_atm"], 99.7))
    field_panel(
        axes[1, 0], xz, arrays["y"] * 1e3, arrays["T"][:, zoom],
        "Temperature", "T (K)", "inferno",
        vmin=np.percentile(arrays["T"], 0.3),
        vmax=np.percentile(arrays["T"], 99.7))
    field_panel(
        axes[1, 1], xz, arrays["y"] * 1e3, arrays["progress"][:, zoom],
        "H2 consumption", "reaction progress", "cividis", 0.0, 1.0)
    fig.suptitle(
        r"Wave-front detail at t = %.3f $\mu$s" % (time * 1e6))
    fig.savefig(outdir / "fields_front_zoom.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 6.2), constrained_layout=True)
    field_panel(
        axes[0, 0], arrays["x"], arrays["y"] * 1e3, arrays["Y_H2"],
        "H2 mass fraction", r"$Y_{H_2}$", "Blues_r")
    field_panel(
        axes[0, 1], arrays["x"], arrays["y"] * 1e3, arrays["Y_O2"],
        "O2 mass fraction", r"$Y_{O_2}$", "Greens_r")
    field_panel(
        axes[1, 0], arrays["x"], arrays["y"] * 1e3, arrays["Y_OH"],
        "OH radical", r"$Y_{OH}$", "plasma")
    field_panel(
        axes[1, 1], arrays["x"], arrays["y"] * 1e3, arrays["Y_H2O"],
        "H2O product", r"$Y_{H_2O}$", "Blues")
    fig.suptitle(r"Species fields at t = %.3f $\mu$s" % (time * 1e6))
    fig.savefig(outdir / "fields_species.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True,
                             constrained_layout=True)
    axes[0].plot(arrays["x"], arrays["p_atm"][mid], label="p")
    axes[0].set_ylabel("p (atm)")
    ax2 = axes[0].twinx()
    ax2.plot(arrays["x"], arrays["T"][mid], color="tab:red", label="T")
    ax2.set_ylabel("T (K)", color="tab:red")
    axes[0].set_title("Centerline pressure and temperature")
    axes[0].grid(alpha=0.25)
    axes[0].set_xlim(zoom_lo, zoom_hi)
    axes[1].plot(arrays["x"], arrays["Y_H2"][mid], label="H2")
    axes[1].plot(arrays["x"], arrays["Y_O2"][mid], label="O2")
    axes[1].plot(arrays["x"], arrays["Y_OH"][mid], label="OH")
    axes[1].plot(arrays["x"], arrays["Y_H2O"][mid], label="H2O")
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("mass fraction")
    axes[1].set_title("Centerline composition")
    axes[1].grid(alpha=0.25)
    axes[1].set_xlim(zoom_lo, zoom_hi)
    axes[1].legend(ncol=4)
    fig.savefig(outdir / "centerline_profiles.png", dpi=220)
    plt.close(fig)

    print("Wrote", outdir)


if __name__ == "__main__":
    main()
