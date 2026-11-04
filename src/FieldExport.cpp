#include "FieldExport.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace amrex;

namespace {

std::string
time_tag(Real time)
{
    std::ostringstream out;
    out << 't' << std::fixed << std::setprecision(12) << time;
    return out.str();
}

std::string
rank_tag(int rank)
{
    std::ostringstream out;
    out << std::setw(6) << std::setfill('0') << rank;
    return out.str();
}

FArrayBox
host_copy(const FArrayBox& source)
{
    const Box& bx = source.box();
    FArrayBox host(bx, source.nComp(), The_Pinned_Arena());
#ifdef AMREX_USE_GPU
    Gpu::copyAsync(
      Gpu::deviceToHost, source.dataPtr(),
      source.dataPtr() + source.size(), host.dataPtr());
    Gpu::streamSynchronize();
#else
    host.copy<RunOn::Host>(source);
#endif
    return host;
}

void
write_vtu(
  const std::filesystem::path& path, const FArrayBox& fab,
  const Box& bx, const Geometry& geom,
  const Vector<std::string>& names, Real time, int precision)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const int nx = hi.x - lo.x + 1;
    const int ny = hi.y - lo.y + 1;
    const int pnx = nx + 1;
    const int pny = ny + 1;
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const long ncells = static_cast<long>(nx) * ny;
    const long npoints = static_cast<long>(pnx) * pny;

    std::ofstream out(path);
    out << std::setprecision(std::max(6, precision));
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
           "byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << npoints
        << "\" NumberOfCells=\"" << ncells << "\">\n"
        << "      <PointData>\n"
        << "      </PointData>\n"
        << "      <CellData Scalars=\"" << names.front() << "\">\n";
    const auto data = fab.const_array();
    for (int n = 0; n < static_cast<int>(names.size()); ++n) {
        out << "        <DataArray type=\"Float64\" Name=\"" << names[n]
            << "\" format=\"ascii\">\n";
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                out << data(lo.x + i, lo.y + j, 0, n) << '\n';
            }
        }
        out << "        </DataArray>\n";
    }
    out << "      </CellData>\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float64\" Name=\"Points\" "
           "NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (int j = 0; j < pny; ++j) {
        for (int i = 0; i < pnx; ++i) {
            const Real x = problo[0] + (lo.x + i) * dx[0];
            const Real y = problo[1] + (lo.y + j) * dx[1];
            out << x << ' ' << y << " 0\n";
        }
    }
    out << "        </DataArray>\n"
        << "      </Points>\n"
        << "      <Cells>\n"
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
           "format=\"ascii\">\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int p00 = j * pnx + i;
            const int p10 = p00 + 1;
            const int p11 = p10 + pnx;
            const int p01 = p00 + pnx;
            out << p00 << ' ' << p10 << ' ' << p11 << ' ' << p01 << '\n';
        }
    }
    out << "        </DataArray>\n"
        << "        <DataArray type=\"Int64\" Name=\"offsets\" "
           "format=\"ascii\">\n";
    for (long c = 0; c < ncells; ++c) {
        out << 4 * (c + 1) << '\n';
    }
    out << "        </DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" "
           "format=\"ascii\">\n";
    for (long c = 0; c < ncells; ++c) {
        out << "9\n";
    }
    out << "        </DataArray>\n"
        << "      </Cells>\n"
        << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";
}

void
write_tecplot_zones(
  std::ostream& out, const FArrayBox& fab, const Box& bx,
  const Geometry& geom, const Vector<std::string>& names,
  Real time, int precision, int lev, int box_index, int rank)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const int nx = hi.x - lo.x + 1;
    const int ny = hi.y - lo.y + 1;
    const int pnx = nx + 1;
    const int pny = ny + 1;
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const int nvars = static_cast<int>(names.size());
    const auto data = fab.const_array();

    out << std::setprecision(std::max(6, precision));
    out << "ZONE T=\"Level " << lev << " box " << box_index
        << " rank " << rank << "\"\n"
        << "  I=" << pnx << " J=" << pny << " K=1\n"
        << "  ZONETYPE=ORDERED DATAPACKING=BLOCK\n"
        << "  VARLOCATION=([3-" << nvars + 2 << "]=CELLCENTERED)\n"
        << "  STRANDID=1 SOLUTIONTIME=" << time << "\n";

    for (int j = 0; j < pny; ++j) {
        for (int i = 0; i < pnx; ++i) {
            out << problo[0] + (lo.x + i) * dx[0] << '\n';
        }
    }
    for (int j = 0; j < pny; ++j) {
        for (int i = 0; i < pnx; ++i) {
            out << problo[1] + (lo.y + j) * dx[1] << '\n';
        }
    }
    for (int n = 0; n < nvars; ++n) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                out << data(lo.x + i, lo.y + j, 0, n) << '\n';
            }
        }
    }
}

void
write_paraview_index(
  const std::filesystem::path& paraview_dir,
  const std::filesystem::path& level_dir, int lev,
  std::ostream& vtm)
{
    Vector<std::string> files;
    if (std::filesystem::exists(level_dir)) {
        for (const auto& entry :
             std::filesystem::directory_iterator(level_dir)) {
            if (entry.is_regular_file() &&
                entry.path().extension() == ".vtu") {
                files.push_back(
                  std::filesystem::relative(
                    entry.path(), paraview_dir).generic_string());
            }
        }
    }
    std::sort(files.begin(), files.end());
    vtm << "    <Block name=\"Level " << lev << "\">\n";
    for (int i = 0; i < static_cast<int>(files.size()); ++i) {
        vtm << "      <DataSet index=\"" << i << "\" file=\""
            << files[i] << "\"/>\n";
    }
    vtm << "    </Block>\n";
}

void
write_time_series_files(
  const std::filesystem::path& paraview_dir)
{
    struct Entry
    {
        Real time;
        std::string name;
    };
    Vector<Entry> entries;
    for (const auto& item :
         std::filesystem::directory_iterator(paraview_dir)) {
        if (!item.is_regular_file() ||
            item.path().extension() != ".vtm") {
            continue;
        }
        const std::string stem = item.path().stem().string();
        if (stem.size() > 1 && stem.front() == 't') {
            try {
                entries.push_back(
                  {static_cast<Real>(std::stod(stem.substr(1))),
                   item.path().filename().string()});
            } catch (const std::exception&) {
            }
        }
    }
    std::sort(
      entries.begin(), entries.end(),
      [](const Entry& a, const Entry& b) { return a.time < b.time; });

    std::ofstream json(paraview_dir / "fields.series");
    json << "{\n"
         << "  \"file-series-version\": \"1.0\",\n"
         << "  \"files\": [\n";
    for (int i = 0; i < static_cast<int>(entries.size()); ++i) {
        json << "    {\"name\": \"" << entries[i].name
             << "\", \"time\": " << entries[i].time << "}";
        json << (i + 1 < static_cast<int>(entries.size()) ? "," : "")
             << '\n';
    }
    json << "  ]\n"
         << "}\n";

    std::ofstream pvd(paraview_dir / "time_series.pvd");
    pvd << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n"
        << "  <Collection>\n";
    for (const auto& entry : entries) {
        pvd << "    <DataSet timestep=\"" << entry.time
            << "\" group=\"\" part=\"0\" file=\"" << entry.name
            << "\"/>\n";
    }
    pvd << "  </Collection>\n"
        << "</VTKFile>\n";
}

} // namespace

namespace detexport {

void
WriteSnapshot(
  const Vector<MultiFab>& fields, const Vector<std::string>& names,
  const Vector<Geometry>& geom, int finest_level, Real time,
  const std::string& output_root, int precision,
  bool write_vtk, bool write_tecplot)
{
    if (!write_vtk && !write_tecplot) {
        return;
    }

    const std::filesystem::path root(output_root);
    const std::filesystem::path paraview_dir = root / "paraview";
    const std::filesystem::path tecplot_dir = root / "tecplot";
    const std::string tag = time_tag(time);
    const std::filesystem::path vtk_snapshot_dir = paraview_dir / tag;
    std::error_code ec;
    if (write_vtk) {
        std::filesystem::create_directories(vtk_snapshot_dir, ec);
    }
    if (write_tecplot) {
        std::filesystem::create_directories(tecplot_dir, ec);
    }
    ParallelDescriptor::Barrier();

    const int rank = ParallelDescriptor::MyProc();
    const int nprocs = ParallelDescriptor::NProcs();
    const std::filesystem::path tec_part =
      tecplot_dir / (tag + ".rank" + rank_tag(rank) + ".part");
    std::ofstream tec;
    if (write_tecplot) {
        tec.open(tec_part);
    }
    tec << std::setprecision(std::max(6, precision));

    for (int lev = 0; lev <= finest_level; ++lev) {
        const std::filesystem::path level_dir =
          vtk_snapshot_dir / ("level_" + std::to_string(lev));
        if (write_vtk) {
            std::filesystem::create_directories(level_dir, ec);
        }
        int box_index = 0;
        for (MFIter mfi(fields[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const FArrayBox host = host_copy(fields[lev][mfi]);
            if (write_vtk) {
                std::ostringstream name;
                name << "level_" << lev << "_rank_" << rank_tag(rank)
                     << "_box_" << std::setw(6) << std::setfill('0')
                     << box_index << ".vtu";
                write_vtu(
                  level_dir / name.str(), host, bx, geom[lev],
                  names, time, precision);
            }
            if (write_tecplot) {
                write_tecplot_zones(
                  tec, host, bx, geom[lev], names,
                  time, precision, lev, box_index, rank);
            }
            ++box_index;
        }
    }
    if (tec.is_open()) {
        tec.close();
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        if (write_vtk) {
            const std::filesystem::path vtm_path =
              paraview_dir / (tag + ".vtm");
            std::ofstream vtm(vtm_path);
            vtm << "<?xml version=\"1.0\"?>\n"
                << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" "
                   "byte_order=\"LittleEndian\">\n"
                << "  <vtkMultiBlockDataSet>\n";
            for (int lev = 0; lev <= finest_level; ++lev) {
                write_paraview_index(
                  paraview_dir,
                  vtk_snapshot_dir / ("level_" + std::to_string(lev)),
                  lev, vtm);
            }
            vtm << "  </vtkMultiBlockDataSet>\n"
                << "</VTKFile>\n";
            write_time_series_files(paraview_dir);
        }

        if (write_tecplot) {
            const std::filesystem::path tec_path =
              tecplot_dir / (tag + ".dat");
            std::ofstream out(tec_path);
            out << std::setprecision(std::max(6, precision))
                << "TITLE = \"AMReX detonation fields at t = "
                << time << " s\"\n"
                << "VARIABLES = \"X\" \"Y\"";
            for (const auto& name : names) {
                out << " \"" << name << "\"";
            }
            out << "\n";
            for (int proc = 0; proc < nprocs; ++proc) {
                const std::filesystem::path part =
                  tecplot_dir /
                  (tag + ".rank" + rank_tag(proc) + ".part");
                if (!std::filesystem::exists(part)) {
                    continue;
                }
                std::ifstream in(part);
                out << in.rdbuf();
                in.close();
                std::filesystem::remove(part, ec);
            }
        }
    }
    ParallelDescriptor::Barrier();
}

} // namespace detexport
