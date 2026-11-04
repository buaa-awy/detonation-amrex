#pragma once

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace detexport {

void WriteSnapshot(
  const amrex::Vector<amrex::MultiFab>& fields,
  const amrex::Vector<std::string>& names,
  const amrex::Vector<amrex::Geometry>& geom,
  int finest_level, amrex::Real time,
  const std::string& output_root, int precision,
  bool write_vtk, bool write_tecplot);

} // namespace detexport
