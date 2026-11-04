#!/usr/bin/env python3
"""Rasterize one solver VTM snapshot without AMReX fextract.

Run this script with ParaView's pvpython:

    pvpython tools/extract_vtm_fields.py snapshot.vtm output.npz

The finest available AMR block is written over coarser blocks.  The resulting
NPZ can be consumed by plot_final_flow.py through its --data option.
"""

import argparse
from pathlib import Path
import re

import numpy as np
import paraview.simple as pvsimple
import paraview.servermanager as servermanager
import vtk
from vtk.util.numpy_support import vtk_to_numpy


FIELD_NAMES = [
    "p", "T", "U_x", "U_y", "rho", "magGradrho", "Qdot", "progress",
    "Y_H2", "Y_H", "Y_O", "Y_O2", "Y_OH", "Y_H2O", "Y_HO2",
    "Y_H2O2", "Y_N2",
]


def visit_grids(dataset, level, callback):
    if isinstance(dataset, vtk.vtkMultiBlockDataSet):
        for index in range(dataset.GetNumberOfBlocks()):
            block = dataset.GetBlock(index)
            if block is not None:
                visit_grids(block, level, callback)
    elif isinstance(dataset, vtk.vtkUnstructuredGrid):
        callback(dataset, level)


def time_from_name(path):
    match = re.search(r"t([0-9.eE+-]+)\.vtm$", path.name)
    if not match:
        raise ValueError(
            "Cannot infer time from VTM name; expected a name like "
            "t0.000070001234.vtm")
    return float(match.group(1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vtm", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--nx", type=int, default=3001)
    parser.add_argument("--ny", type=int, default=401)
    args = parser.parse_args()

    # Do not call resolve(): on Windows it may expand an ASCII junction back
    # to a non-ASCII target path, which ParaView's Python file reader rejects.
    source = pvsimple.OpenDataFile(str(args.vtm.absolute()))
    source.UpdatePipeline()
    dataset = servermanager.Fetch(source)
    if not isinstance(dataset, vtk.vtkMultiBlockDataSet):
        raise RuntimeError("Expected a VTK multi-block data set")

    grids = []
    for level in range(dataset.GetNumberOfBlocks()):
        block = dataset.GetBlock(level)
        if block is not None:
            visit_grids(block, level, lambda grid, lev: grids.append((grid, lev)))
    if not grids:
        raise RuntimeError("VTM snapshot has no unstructured-grid blocks")
    xmin = min(grid.GetBounds()[0] for grid, _ in grids)
    xmax = max(grid.GetBounds()[1] for grid, _ in grids)
    ymin = min(grid.GetBounds()[2] for grid, _ in grids)
    ymax = max(grid.GetBounds()[3] for grid, _ in grids)
    nx = max(2, args.nx)
    ny = max(2, args.ny)
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    values = {
        name: np.full((ny, nx), np.nan, dtype=np.float64)
        for name in FIELD_NAMES
    }
    levels = np.full((ny, nx), -1, dtype=np.int16)

    def write_grid(grid, level):
        centers = vtk.vtkCellCenters()
        centers.SetInputData(grid)
        centers.Update()
        sampled = centers.GetOutput()
        points = vtk_to_numpy(sampled.GetPoints().GetData())
        cx = np.unique(np.round(points[:, 0], 12))
        cy = np.unique(np.round(points[:, 1], 12))
        if cx.size * cy.size != points.shape[0]:
            raise RuntimeError("VTM block is not a structured 2-D grid")
        bounds = grid.GetBounds()
        tx0 = np.searchsorted(x, bounds[0], side="left")
        tx1 = np.searchsorted(x, bounds[1], side="right")
        ty0 = np.searchsorted(y, bounds[2], side="left")
        ty1 = np.searchsorted(y, bounds[3], side="right")
        ix = np.abs(x[tx0:tx1, None] - cx[None, :]).argmin(axis=1)
        iy = np.abs(y[ty0:ty1, None] - cy[None, :]).argmin(axis=1)
        array_data = sampled.GetPointData()
        for name in FIELD_NAMES:
            array = array_data.GetArray(name)
            if array is None:
                raise RuntimeError("VTM snapshot is missing field " + name)
            block_values = vtk_to_numpy(array).reshape(cy.size, cx.size)
            values[name][ty0:ty1, tx0:tx1] = (
                block_values[iy[:, None], ix[None, :]])
        levels[ty0:ty1, tx0:tx1] = level

    # The root VTM block order is AMR level order.  Later, finer levels
    # overwrite coarse values in the raster.
    for grid, level in grids:
        write_grid(grid, level)

    for name, field in values.items():
        if np.isnan(field).any():
            missing = np.argwhere(np.isnan(field))
            preview = ", ".join(
                "(%.6g, %.6g)" % (x[i], y[j])
                for j, i in missing[:8])
            raise RuntimeError(
                "Rasterization left %d uncovered cells in %s: %s"
                % (missing.shape[0], name, preview))
    values["p_atm"] = values["p"] / 101325.0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        x=x,
        y=y,
        time_us=np.array(time_from_name(args.vtm) * 1.0e6),
        level=levels,
        **values
    )
    print("Wrote", args.output)


if __name__ == "__main__":
    main()
