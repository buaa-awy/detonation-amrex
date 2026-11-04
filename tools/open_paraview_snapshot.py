"""Open one AMReX VTK snapshot in the ParaView desktop application."""

import os
import sys

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    OpenDataFile,
    Render,
    ResetCamera,
    SetActiveSource,
    Show,
)


path = os.environ.get("DETONATION_VTM")
if not path and len(sys.argv) > 1:
    path = sys.argv[1]
if not path:
    raise RuntimeError(
        "Set DETONATION_VTM or pass the .vtm path to the ParaView script.")
source = OpenDataFile(path)
source.UpdatePipeline()
SetActiveSource(source)

view = GetActiveViewOrCreate("RenderView")
display = Show(source, view)
display.Representation = "Surface"
ColorBy(display, ("CELLS", "T"))
display.SetScalarBarVisibility(view, True)
ResetCamera(view)
Render(view)
