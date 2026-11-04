#include "Detonation.H"
#include "FieldExport.H"

#include <AMReX_BC_TYPES.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace amrex;
using namespace amrex::literals;
using namespace det;

namespace {

struct DetFill
{
    AMREX_GPU_DEVICE
    void
    operator()(
      const IntVect&, Array4<Real> const&, int, int,
      GeometryData const&, Real, const BCRec*, int, int) const
    {
    }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
minmod(Real a, Real b)
{
    if (a * b <= 0.0_rt) {
        return 0.0_rt;
    }
    return (std::abs(a) < std::abs(b)) ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_reconstructed(
  const Array4<Real const>& q, const Array4<Real const>& sl,
  int i, int j, int dir, Real sign, Real* out)
{
    bool bad = false;
    for (int n = 0; n < det::NQ; ++n) {
        const Real slope = sl(i, j, 0, n + dir * det::NQ);
        out[n] = q(i, j, 0, n) + sign * 0.5_rt * slope;
    }
    bad = out[det::QRHO] <= 0.0_rt || out[det::QPRES] <= 0.0_rt;
    if (bad) {
        for (int n = 0; n < det::NQ; ++n) {
            out[n] = q(i, j, 0, n);
        }
    }
    out[det::QTEMP] = std::max(out[det::QTEMP], 250.0_rt);
    out[det::QCS] = std::max(out[det::QCS], 1.0_rt);
    out[det::QGAM] = std::max(out[det::QGAM], 1.000001_rt);
    Real sum = 0.0_rt;
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        out[det::QY + k] = std::max(out[det::QY + k], 0.0_rt);
        sum += out[det::QY + k];
    }
    sum = std::max(sum, 1.0e-30_rt);
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        out[det::QY + k] /= sum;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real
total_energy_from_primitive(
  const Real* q, const Real* Y)
{
    Real e = 0.0_rt;
    Real cp = 0.0_rt;
    Real Rmix = 0.0_rt;
    mixture_thermo(Y, q[det::QTEMP], e, cp, Rmix);
    return q[det::QRHO] *
           (e + 0.5_rt * (q[det::QU] * q[det::QU] +
                           q[det::QV] * q[det::QV]));
}

void
fill_current_impl(
  const Geometry& cgeom, const Geometry& fgeom,
  const Vector<BCRec>& bcs,
  const Vector<MultiFab>& source, int lev, Real time,
  MultiFab& data, int ngrow, const IntVect& rr)
{
    Vector<MultiFab*> smf;
    Vector<Real> stime;
    if (lev == 0) {
        smf.push_back(const_cast<MultiFab*>(&source[0]));
        stime.push_back(time);
    } else {
        smf.push_back(const_cast<MultiFab*>(&source[lev - 1]));
        smf.push_back(const_cast<MultiFab*>(&source[lev]));
        stime.push_back(time);
        stime.push_back(time);
    }

    if (Gpu::inLaunchRegion()) {
        GpuBndryFuncFab<DetFill> bndry_func(DetFill{});
        if (lev == 0) {
            PhysBCFunct<GpuBndryFuncFab<DetFill>> physbc(
              fgeom, bcs, bndry_func);
            FillPatchSingleLevel(
              data, time, smf, stime, 0, 0, det::NCOMP, fgeom, physbc, 0);
        } else {
            PhysBCFunct<GpuBndryFuncFab<DetFill>> cphysbc(
              cgeom, bcs, bndry_func);
            PhysBCFunct<GpuBndryFuncFab<DetFill>> fphysbc(
              fgeom, bcs, bndry_func);
            FillPatchTwoLevels(
              data, time, {smf[0]}, {stime[0]}, {smf[1]}, {stime[1]},
              0, 0, det::NCOMP, cgeom, fgeom, cphysbc, 0, fphysbc, 0,
              rr, &cell_cons_interp, bcs, 0);
        }
    } else {
        CpuBndryFuncFab bndry_func(nullptr);
        if (lev == 0) {
            PhysBCFunct<CpuBndryFuncFab> physbc(fgeom, bcs, bndry_func);
            FillPatchSingleLevel(
              data, time, smf, stime, 0, 0, det::NCOMP, fgeom, physbc, 0);
        } else {
            PhysBCFunct<CpuBndryFuncFab> cphysbc(cgeom, bcs, bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(fgeom, bcs, bndry_func);
            FillPatchTwoLevels(
              data, time, {smf[0]}, {stime[0]}, {smf[1]}, {stime[1]},
              0, 0, det::NCOMP, cgeom, fgeom, cphysbc, 0, fphysbc, 0,
              rr, &cell_cons_interp, bcs, 0);
        }
    }
}

void
average_faces(
  Vector<Array<MultiFab, AMREX_SPACEDIM>>& fluxes,
  const Vector<Geometry>& geom, const Vector<IntVect>& rr, int finest)
{
    for (int lev = finest; lev > 0; --lev) {
        average_down_faces(
          GetArrOfConstPtrs(fluxes[lev]), GetArrOfPtrs(fluxes[lev - 1]),
          rr[lev - 1], geom[lev - 1]);
    }
}

} // namespace

Detonation::Detonation()
{
    ReadParameters();
    const int nlevels = max_level + 1;
    state.resize(nlevels);
    chem_temperature.resize(nlevels);
    chem_hnext.resize(nlevels);
    chem_activity.resize(nlevels);
    fine_mask.resize(nlevels);
    t_new.resize(nlevels, 0.0_rt);
    t_old.resize(nlevels, -1.0e30_rt);
    dt.resize(nlevels, 1.0e30_rt);
    istep.resize(nlevels, 0);
    fillpatcher.resize(nlevels);
    flux_register.resize(nlevels);

    bcs.resize(det::NCOMP);
    for (int n = 0; n < det::NCOMP; ++n) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bcs[n].setLo(dir, BCType::foextrap);
            bcs[n].setHi(dir, BCType::foextrap);
        }
        bcs[n].setLo(0, (n == det::IMX) ? BCType::reflect_odd
                                        : BCType::reflect_even);
        bcs[n].setHi(0, BCType::foextrap);
        if (AMREX_SPACEDIM > 1) {
            bcs[n].setLo(1, BCType::int_dir);
            bcs[n].setHi(1, BCType::int_dir);
        }
    }
    wall_start = ParallelDescriptor::second();
}

void Detonation::ReadParameters()
{
    ParmParse pp;
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);

    ParmParse amr_pp("amr");
    amr_pp.query("regrid_int", regrid_int);
    amr_pp.query("plot_int", plot_int);
    amr_pp.query("chk_int", chk_int);
    amr_pp.query("plot_time_int", plot_time_int);
    amr_pp.query("chk_time_int", chk_time_int);
    amr_pp.query("plot_file", plot_file);
    amr_pp.query("chk_file", chk_file);
    amr_pp.query("restart", restart_chkfile);
    amr_pp.query("do_reflux", do_reflux);
    amr_pp.query("verbose", verbose);
    amr_pp.query("plot_vtk", export_vtk);
    amr_pp.query("plot_tecplot", export_tecplot);
    amr_pp.query("plot_export_precision", export_precision);

    ParmParse det_pp("det");
    det_pp.query("cfl", cfl);
    det_pp.query("change_max", change_max);
    det_pp.query("x_igniter", x_igniter);
    det_pp.query("T_igniter", T_igniter);
    det_pp.query("p_igniter", p_igniter);
    det_pp.query("T_fuel", T_fuel);
    det_pp.query("p_fuel", p_fuel);
    det_pp.query("perturb_amp", perturb_amp);
    det_pp.query("perturb_waves", perturb_waves);
    det_pp.query("Y_H2", Y_H2);
    det_pp.query("Y_O2", Y_O2);
    det_pp.query("grad_p", grad_p);
    det_pp.query("grad_T", grad_T);
    det_pp.query("grad_Y", grad_Y);
    det_pp.query("chem_activity", chem_activity_tol);
    det_pp.query("front_pre", front_pre);
    det_pp.query("front_post", front_post);
    det_pp.query("refine_zone_lo", refine_zone_lo);
    det_pp.query("refine_zone_hi", refine_zone_hi);
    det_pp.query("refine_zone_until", refine_zone_until);
    det_pp.query("front_regrid_distance", front_regrid_distance);
    det_pp.query("chem_Tmin", Tmin_chem);
    det_pp.query("chem_hmax", hmax_chem);
    det_pp.query("chem_tol_y", tol_y_chem);
    det_pp.query("chem_tol_T", tol_T_chem);
    det_pp.query("progress_every", progress_int);
    det_pp.query("status_every", status_int);
    det_pp.query("chem_stride", chem_stride);
    det_pp.query("status_file", status_file);
    det_pp.query("status_text_file", status_text_file);
    det_pp.query("allow_legacy_restart", allow_legacy_restart);
    front_x_cached = x_igniter;
    next_plot_time = plot_time_int;
    next_chk_time = chk_time_int;
}

void Detonation::Evolve()
{
    Real time = t_new[0];
    while (istep[0] < max_step && time < stop_time) {
        const Real finest_dx = Geom(finest_level).CellSize(0);
        const bool regrid_by_step =
          (regrid_int > 0 && istep[0] % regrid_int == 0);
        const bool regrid_by_front =
          (last_regrid_front >= 0.0_rt &&
           std::abs(front_x_cached - last_regrid_front) >=
             std::max(front_regrid_distance, 4.0_rt * finest_dx));
        if (max_level > 0 && (regrid_by_step || regrid_by_front)) {
            if (verbose) {
                Print() << "REGRID_SYNC step=" << istep[0] << '\n';
            }
            // Level storage is rebuilt during regrid.  Complete all GPU work
            // from the preceding step before releasing/reallocating it.
            Gpu::streamSynchronize();
            if (verbose) {
                Print() << "REGRID_BEGIN step=" << istep[0] << '\n';
            }
            regrid(0, time);
            last_regrid_front = front_x_cached;
            UpdateFineMasks();
            if (verbose) {
                Print() << "REGRID_END step=" << istep[0] << '\n';
            }
        }
        AverageDown();
        const Real dt0 = ComputeDt();

        const int chem_stride_loc = std::max(1, chem_stride);
        const bool begin_block = (istep[0] % chem_stride_loc == 0);
        if (begin_block) {
            chem_block_dt = 0.0_rt;
            ChemistryHalfStep(
              0.5_rt * chem_stride_loc * dt0);
        }
        HydroRK2(time, dt0);
        chem_block_dt += dt0;
        if ((istep[0] + 1) % chem_stride_loc == 0) {
            ChemistryHalfStep(0.5_rt * chem_block_dt);
        }
        AverageDown();
        if (!CheckFiniteState("step", time + dt0)) {
            Abort("Non-finite conserved state");
        }

        time += dt0;
        for (int lev = 0; lev <= finest_level; ++lev) {
            t_old[lev] = t_new[lev];
            t_new[lev] = time;
            dt[lev] = dt0;
            ++istep[lev];
        }

        const Real wall = ParallelDescriptor::second() - wall_start;
        const Real frac = (stop_time > 0.0_rt)
                            ? std::min(1.0_rt, time / stop_time) : 1.0_rt;
        const Real eta = EstimateEta(time, wall, false);
        if (progress_int > 0 && istep[0] % progress_int == 0) {
            Print() << "PROGRESS step=" << istep[0]
                    << " t_us=" << time * 1.0e6
                    << " target_us=" << stop_time * 1.0e6
                    << " percent=" << frac * 100.0_rt
                    << " dt_s=" << dt0
                    << " finest_level=" << finest_level
                    << " wall_s=" << wall
                    << " eta_s=" << eta << '\n';
            Print() << "CONSERVED step=" << istep[0]
                    << " mass=" << state[0].sum(det::IRHO)
                    << " energy=" << state[0].sum(det::IENE)
                    << '\n';
        }
        if (status_int > 0 && istep[0] % status_int == 0) {
            WriteStatusFile(time, dt0, wall, eta, false);
        }
        const bool plot_by_step =
          (plot_int > 0 && istep[0] % plot_int == 0);
        const bool plot_by_time =
          (plot_time_int > 0.0_rt && time >= next_plot_time);
        if (plot_by_step || plot_by_time) {
            WritePlotFile();
            last_plot_time = time;
            if (plot_by_time) {
                do {
                    next_plot_time += plot_time_int;
                } while (next_plot_time <= time);
            }
        }
        const bool chk_by_step =
          (chk_int > 0 && istep[0] % chk_int == 0);
        const bool chk_by_time =
          (chk_time_int > 0.0_rt && time >= next_chk_time);
        if (chk_by_step || chk_by_time) {
            WriteCheckpointFile();
            last_chk_time = time;
            if (chk_by_time) {
                do {
                    next_chk_time += chk_time_int;
                } while (next_chk_time <= time);
            }
        }
    }
    const Real output_eps = 1.0e-12_rt * std::max(stop_time, 1.0_rt);
    if ((plot_int > 0 || plot_time_int > 0.0_rt) &&
        last_plot_time < time - output_eps) {
        WritePlotFile();
        last_plot_time = time;
    }
    if ((chk_int > 0 || chk_time_int > 0.0_rt) &&
        last_chk_time < time - output_eps) {
        WriteCheckpointFile();
        last_chk_time = time;
    }
    const Real wall = ParallelDescriptor::second() - wall_start;
    WriteStatusFile(time, dt[0], wall, 0.0_rt, true);
}

Real Detonation::EstimateEta(Real time, Real wall, bool finished)
{
    if (finished || time <= 0.0_rt) {
        return 0.0_rt;
    }
    Real estimate = -1.0_rt;
    if (eta_last_time >= 0.0_rt && time > eta_last_time &&
        wall > eta_last_wall) {
        const Real recent_rate =
          (wall - eta_last_wall) / (time - eta_last_time);
        estimate = recent_rate * std::max(0.0_rt, stop_time - time);
        eta_smoothed =
          (eta_smoothed >= 0.0_rt)
            ? 0.75_rt * eta_smoothed + 0.25_rt * estimate
            : estimate;
    } else {
        estimate = wall * std::max(0.0_rt, stop_time - time) / time;
        eta_smoothed = estimate;
    }
    eta_last_time = time;
    eta_last_wall = wall;
    return eta_smoothed;
}

void Detonation::InitData()
{
    if (restart_chkfile.empty()) {
        InitFromScratch(0.0_rt);
        AverageDown();
        WritePlotFile();
        last_plot_time = 0.0_rt;
        if (chk_int > 0 || chk_time_int > 0.0_rt) {
            WriteCheckpointFile();
            last_chk_time = 0.0_rt;
        }
    } else {
        ReadCheckpointFile();
        // Rebuild AMR boxes under the current tagging/max_grid_size rules.
        // Older checkpoints may contain boxes that are not slab-aligned in the
        // periodic transverse direction.
        regrid(0, t_new[0]);
        AverageDown();
        WritePlotFile();
        last_plot_time = t_new[0];
        last_chk_time = t_new[0];
        while (plot_time_int > 0.0_rt &&
               next_plot_time <= t_new[0]) {
            next_plot_time += plot_time_int;
        }
        while (chk_time_int > 0.0_rt &&
               next_chk_time <= t_new[0]) {
            next_chk_time += chk_time_int;
        }
    }
    UpdateFineMasks();
    Print() << "INITIAL_CONSERVED mass=" << state[0].sum(det::IRHO)
            << " energy=" << state[0].sum(det::IENE) << '\n';
    const Real wall = ParallelDescriptor::second() - wall_start;
    WriteStatusFile(t_new[0], 0.0_rt, wall, -1.0_rt, false);
    last_regrid_front = front_x_cached;
    front_x_previous = front_x_cached;
    front_time_previous = t_new[0];
}

void Detonation::InitLevelData(int lev)
{
    MultiFab& U = state[lev];
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();
    const Real Ly = Geom(lev).ProbHi(1) - Geom(lev).ProbLo(1);
    const Real x_ign = x_igniter;
    const Real T_ign = T_igniter;
    const Real p_ign = p_igniter;
    const Real T_fuel_loc = T_fuel;
    const Real p_fuel_loc = p_fuel;
    const Real amp = perturb_amp;
    const int waves = perturb_waves;
    const Real y_h2 = Y_H2;
    const Real y_o2 = Y_O2;

    for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u = U.array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real x = problo[0] + (i + 0.5_rt) * dx[0];
              const Real y = problo[1] + (j + 0.5_rt) * dx[1];
              Real Y[detmech::NUM_SPECIES] = {0.0_rt};
              Y[0] = y_h2;
              Y[3] = y_o2;
              Real T = (x <= x_ign) ? T_ign : T_fuel_loc;
              const Real p = (x <= x_ign) ? p_ign : p_fuel_loc;
              if (x <= x_ign && amp != 0.0_rt) {
                  const Real phase =
                    2.0_rt * 3.14159265358979323846_rt *
                    waves * y / Ly;
                  T *= 1.0_rt + amp * std::sin(phase);
              }
              Real e = 0.0_rt;
              Real cp = 0.0_rt;
              Real Rmix = 0.0_rt;
              mixture_thermo(Y, T, e, cp, Rmix);
              const Real rho = p / (Rmix * T);
              u(i, j, k, det::IRHO) = rho;
              u(i, j, k, det::IMX) = 0.0_rt;
              u(i, j, k, det::IMY) = 0.0_rt;
              u(i, j, k, det::IENE) = rho * e;
              for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                  u(i, j, k, det::IFY + s) = rho * Y[s];
              }
          });
    }
}

void Detonation::InitializeTemperatureGuess(int lev)
{
    chem_temperature[lev].setVal(T_fuel);
    chem_hnext[lev].setVal(hmax_chem);
    chem_activity[lev].setVal(0.0_rt);
    for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u = state[lev].const_array(mfi);
        auto ta = chem_temperature[lev].array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real Uloc[det::NCOMP];
              Real Y[detmech::NUM_SPECIES];
              for (int n = 0; n < det::NCOMP; ++n) {
                  Uloc[n] = u(i, j, k, n);
              }
              Real Tguess = 1000.0_rt;
              Real rho, vx, vy, p, T, cs, gamma;
              primitive_from_conserved(
                Uloc, Tguess, rho, vx, vy, p, T, cs, gamma, Y);
              ta(i, j, k) = T;
          });
    }
}

void Detonation::UpdateFineMasks()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (lev < finest_level) {
            fine_mask[lev] = makeFineMask(
              state[lev], state[lev + 1], IntVect(0),
              refRatio(lev), geom[lev].periodicity(), 0, 1);
        } else {
            fine_mask[lev].define(
              grids[lev], dmap[lev], 1, 0);
            fine_mask[lev].setVal(0);
        }
    }
}

void Detonation::MakeNewLevelFromScratch(
  int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    state[lev].define(ba, dm, det::NCOMP, 3);
    chem_temperature[lev].define(ba, dm, 1, 0);
    chem_hnext[lev].define(ba, dm, 1, 0);
    chem_activity[lev].define(ba, dm, 1, 0);
    t_new[lev] = time;
    t_old[lev] = time - 1.0e30_rt;
    dt[lev] = 1.0e30_rt;
    InitLevelData(lev);
    InitializeTemperatureGuess(lev);
}

void Detonation::MakeNewLevelFromCoarse(
  int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    state[lev].define(ba, dm, det::NCOMP, 3);
    chem_temperature[lev].define(ba, dm, 1, 0);
    chem_hnext[lev].define(ba, dm, 1, 0);
    chem_activity[lev].define(ba, dm, 1, 0);
    if (do_reflux && lev > 0) {
        flux_register[lev] = std::make_unique<FluxRegister>(
          ba, dm, refRatio(lev - 1), lev, det::NCOMP);
    }
    t_new[lev] = time;
    t_old[lev] = time - 1.0e30_rt;
    dt[lev] = 1.0e30_rt;
    FillCoarsePatch(lev, time, state[lev]);
    InitializeTemperatureGuess(lev);
}

void Detonation::RemakeLevel(
  int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    MultiFab new_state(ba, dm, det::NCOMP, 3);
    MultiFab new_temperature(ba, dm, 1, 0);
    MultiFab new_hnext(ba, dm, 1, 0);
    MultiFab new_activity(ba, dm, 1, 0);
    FillCurrent(lev, time, state, new_state, 3);
    Gpu::streamSynchronize();
    state[lev] = std::move(new_state);
    chem_temperature[lev] = std::move(new_temperature);
    chem_hnext[lev] = std::move(new_hnext);
    chem_activity[lev] = std::move(new_activity);
    if (do_reflux && lev > 0) {
        flux_register[lev] = std::make_unique<FluxRegister>(
          ba, dm, refRatio(lev - 1), lev, det::NCOMP);
    }
    t_new[lev] = time;
    dt[lev] = std::max(dt[lev], 0.0_rt);
    InitializeTemperatureGuess(lev);
}

void Detonation::ClearLevel(int lev)
{
    state[lev].clear();
    chem_temperature[lev].clear();
    chem_hnext[lev].clear();
    chem_activity[lev].clear();
    fine_mask[lev].clear();
    flux_register[lev].reset();
}

void Detonation::FillCoarsePatch(int lev, Real time, MultiFab& data)
{
    Vector<MultiFab*> cmf{&state[lev - 1]};
    Vector<Real> ctime{time};
    if (Gpu::inLaunchRegion()) {
        GpuBndryFuncFab<DetFill> func(DetFill{});
        PhysBCFunct<GpuBndryFuncFab<DetFill>> cbc(geom[lev - 1], bcs, func);
        PhysBCFunct<GpuBndryFuncFab<DetFill>> fbc(geom[lev], bcs, func);
        InterpFromCoarseLevel(
          data, time, *cmf[0], 0, 0, det::NCOMP, geom[lev - 1], geom[lev],
          cbc, 0, fbc, 0, refRatio(lev - 1), &cell_cons_interp, bcs, 0);
    } else {
        CpuBndryFuncFab func(nullptr);
        PhysBCFunct<CpuBndryFuncFab> cbc(geom[lev - 1], bcs, func);
        PhysBCFunct<CpuBndryFuncFab> fbc(geom[lev], bcs, func);
        InterpFromCoarseLevel(
          data, time, *cmf[0], 0, 0, det::NCOMP, geom[lev - 1], geom[lev],
          cbc, 0, fbc, 0, refRatio(lev - 1), &cell_cons_interp, bcs, 0);
    }
}

void Detonation::FillCurrent(
  int lev, Real time, const Vector<MultiFab>& source,
  MultiFab& data, int ngrow)
{
    if (lev == 0) {
        fill_current_impl(
          geom[lev], geom[lev], bcs, source, lev, time, data, ngrow,
          IntVect(1));
    } else {
        fill_current_impl(
          geom[lev - 1], geom[lev], bcs, source, lev, time, data, ngrow,
          refRatio(lev - 1));
    }
    if (geom[lev].isAnyPeriodic()) {
        Gpu::streamSynchronize();
        data.FillBoundary(
          0, det::NCOMP, IntVect(ngrow),
          geom[lev].periodicity(), false);
        Gpu::streamSynchronize();
    }
}

void Detonation::AverageDown()
{
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        average_down(
          state[lev + 1], state[lev], geom[lev + 1], geom[lev],
          0, det::NCOMP, refRatio(lev));
        if (chem_temperature[lev + 1].ok() &&
            chem_temperature[lev].ok()) {
            average_down(
              chem_temperature[lev + 1], chem_temperature[lev],
              geom[lev + 1], geom[lev], 0, 1, refRatio(lev));
            average_down(
              chem_hnext[lev + 1], chem_hnext[lev],
              geom[lev + 1], geom[lev], 0, 1, refRatio(lev));
        }
    }
}

Real Detonation::EstTimeStep(int lev)
{
    MultiFab spectral(grids[lev], dmap[lev], 1, 0);
    const auto dx = geom[lev].CellSizeArray();
    for (MFIter mfi(spectral, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto sf = spectral.array(mfi);
        auto u = state[lev].const_array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real rho = std::max(u(i, j, k, det::IRHO), 1.0e-14_rt);
              const Real vx = u(i, j, k, det::IMX) / rho;
              const Real vy = u(i, j, k, det::IMY) / rho;
              Real Y[detmech::NUM_SPECIES];
              Real Tguess = 1000.0_rt;
              Real rr = rho;
              Real uu = 0.0_rt;
              Real vv = 0.0_rt;
              Real p = 0.0_rt;
              Real T = 0.0_rt;
              Real cs = 0.0_rt;
              Real gamma = 1.4_rt;
              Real Uloc[det::NCOMP];
              for (int n = 0; n < det::NCOMP; ++n) {
                  Uloc[n] = u(i, j, k, n);
              }
              primitive_from_conserved(
                Uloc, Tguess, rr, uu, vv, p, T, cs, gamma, Y);
              sf(i, j, k) =
                (std::abs(vx) + cs) / dx[0] +
                (std::abs(vy) + cs) / dx[1];
          });
    }
    const Real smax = spectral.norminf(0, 0, true);
    return (smax > 0.0_rt) ? cfl / smax
                           : std::numeric_limits<Real>::max();
}

Real Detonation::ComputeDt()
{
    Real dt0 = std::numeric_limits<Real>::max();
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt0 = std::min(dt0, EstTimeStep(lev));
    }
    ParallelDescriptor::ReduceRealMin(dt0);
    if (t_new[0] > 0.0_rt && dt[0] < 1.0e29_rt) {
        dt0 = std::min(dt0, change_max * dt[0]);
    }
    const Real eps = 1.0e-12_rt * std::max(stop_time, 1.0_rt);
    dt0 = std::min(dt0, stop_time - t_new[0] + eps);
    return std::max(dt0, 1.0e-20_rt);
}

void Detonation::ChemistryHalfStep(Real dt_half)
{
    if (dt_half <= 0.0_rt) {
        return;
    }
    const Real Tmin = Tmin_chem;
    const Real hmax = hmax_chem;
    const Real tol_y = tol_y_chem;
    const Real tol_T = tol_T_chem;
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto u = state[lev].array(mfi);
            auto ta = chem_temperature[lev].array(mfi);
            auto ha = chem_hnext[lev].array(mfi);
            auto aa = chem_activity[lev].array(mfi);
            auto fm = fine_mask[lev].const_array(mfi);
            ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  if (fm(i, j, k) != 0) {
                      aa(i, j, k) = 0.0_rt;
                      return;
                  }
                  const Real rho = std::max(u(i, j, k, det::IRHO), 1.0e-12_rt);
                  Real Y[detmech::NUM_SPECIES];
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      Y[s] = std::max(u(i, j, k, det::IFY + s), 0.0_rt) / rho;
                  }
                  Real sr = 0.0_rt;
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      sr += Y[s];
                  }
                  sr = std::max(sr, 1.0e-30_rt);
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      Y[s] /= sr;
                  }
                  const Real vx = u(i, j, k, det::IMX) / rho;
                  const Real vy = u(i, j, k, det::IMY) / rho;
                  const Real e0 =
                    (u(i, j, k, det::IENE) -
                     0.5_rt * rho * (vx * vx + vy * vy)) / rho;
                  Real Tguess = ta(i, j, k);
                  if (!std::isfinite(Tguess) ||
                      Tguess < 250.0_rt || Tguess > 6500.0_rt) {
                      Tguess = 1000.0_rt;
                  }
                  Real T = temperature_from_energy(e0, Y, Tguess);
                  Real hnext = ha(i, j, k);
                  Real activity = 0.0_rt;
                  chemistry_step(
                    rho, e0, dt_half, Tmin, hmax, tol_y, tol_T, Y, T,
                    hnext, activity);
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      u(i, j, k, det::IFY + s) = rho * Y[s];
                  }
                  ta(i, j, k) = T;
                  ha(i, j, k) = hnext;
                  aa(i, j, k) = activity;
              });
        }
    }
}

void Detonation::ComputeFluxes(
  int lev, Real time, const Vector<MultiFab>& source,
  Array<MultiFab, AMREX_SPACEDIM>& fluxes)
{
    MultiFab sborder(grids[lev], dmap[lev], det::NCOMP, 3);
    MultiFab tborder(grids[lev], dmap[lev], 1, 3);
    FillCurrent(lev, time, source, sborder, 3);
    tborder.setVal(1000.0_rt);
    MultiFab::Copy(tborder, chem_temperature[lev], 0, 0, 1, 0);
    if (geom[lev].isAnyPeriodic()) {
        tborder.FillBoundary(geom[lev].periodicity());
    }

    for (MFIter mfi(sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Box& qbox = amrex::grow(bx, 2);
        // The WSL CUDA driver is not reliable with repeated stream-ordered
        // frees while AMR level storage is rebuilt. The regular arena keeps
        // these short-lived work arrays until their GPU work has completed.
        FArrayBox qtmp(qbox, det::NQ, The_Arena());
        FArrayBox sltmp(amrex::grow(bx, 1), 2 * det::NQ, The_Arena());
        Elixir e1 = qtmp.elixir();
        Elixir e2 = sltmp.elixir();
        auto q = qtmp.array();
        auto sl = sltmp.array();
        auto u = sborder.const_array(mfi);
        auto ta = tborder.const_array(mfi);

        ParallelFor(
          qbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real Uloc[det::NCOMP];
              for (int n = 0; n < det::NCOMP; ++n) {
                  Uloc[n] = u(i, j, k, n);
              }
              Real Y[detmech::NUM_SPECIES];
              Real Tguess = ta(i, j, k);
              if (!std::isfinite(Tguess) ||
                  Tguess < 250.0_rt || Tguess > 6500.0_rt) {
                  Tguess = 1000.0_rt;
              }
              Real rho, vx, vy, p, T, cs, gamma;
              primitive_from_conserved(
                Uloc, Tguess, rho, vx, vy, p, T, cs, gamma, Y);
              q(i, j, k, det::QRHO) = rho;
              q(i, j, k, det::QU) = vx;
              q(i, j, k, det::QV) = vy;
              q(i, j, k, det::QPRES) = p;
              q(i, j, k, det::QTEMP) = T;
              q(i, j, k, det::QCS) = cs;
              q(i, j, k, det::QGAM) = gamma;
              for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                  q(i, j, k, det::QY + s) = Y[s];
              }
          });
        ParallelFor(
          amrex::grow(bx, 1),
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < det::NQ; ++n) {
                  sl(i, j, k, n) =
                    minmod(q(i, j, k, n) - q(i - 1, j, k, n),
                           q(i + 1, j, k, n) - q(i, j, k, n));
                  sl(i, j, k, det::NQ + n) =
                    minmod(q(i, j, k, n) - q(i, j - 1, k, n),
                           q(i, j + 1, k, n) - q(i, j, k, n));
              }
          });

        auto fx = fluxes[0].array(mfi);
        auto fy = fluxes[1].array(mfi);
        ParallelFor(
          mfi.nodaltilebox(0),
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real qL[det::NQ], qR[det::NQ];
              fill_reconstructed(q, sl, i - 1, j, 0, 1.0_rt, qL);
              fill_reconstructed(q, sl, i, j, 0, -1.0_rt, qR);
              Real YL[detmech::NUM_SPECIES], YR[detmech::NUM_SPECIES];
              for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                  YL[s] = qL[det::QY + s];
                  YR[s] = qR[det::QY + s];
              }
              const Real EL = total_energy_from_primitive(qL, YL);
              const Real ER = total_energy_from_primitive(qR, YR);
              Real F[det::NCOMP];
              hllc_flux(
                0, qL[det::QCS], qR[det::QCS],
                qL[det::QRHO], qL[det::QU], qL[det::QV],
                qL[det::QPRES], EL, YL,
                qR[det::QRHO], qR[det::QU], qR[det::QV],
                qR[det::QPRES], ER, YR, F);
              for (int n = 0; n < det::NCOMP; ++n) {
                  fx(i, j, k, n) = F[n];
              }
          });

        ParallelFor(
          mfi.nodaltilebox(1),
          [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real qL[det::NQ], qR[det::NQ];
              fill_reconstructed(q, sl, i, j - 1, 1, 1.0_rt, qL);
              fill_reconstructed(q, sl, i, j, 1, -1.0_rt, qR);
              Real YL[detmech::NUM_SPECIES], YR[detmech::NUM_SPECIES];
              for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                  YL[s] = qL[det::QY + s];
                  YR[s] = qR[det::QY + s];
              }
              const Real EL = total_energy_from_primitive(qL, YL);
              const Real ER = total_energy_from_primitive(qR, YR);
              Real F[det::NCOMP];
              hllc_flux(
                1, qL[det::QCS], qR[det::QCS],
                qL[det::QRHO], qL[det::QU], qL[det::QV],
                qL[det::QPRES], EL, YL,
                qR[det::QRHO], qR[det::QU], qR[det::QV],
                qR[det::QPRES], ER, YR, F);
              for (int n = 0; n < det::NCOMP; ++n) {
                  fy(i, j, k, n) = F[n];
              }
          });
    }
}

void Detonation::UpdateWithFluxes(
  int lev, Real dt_lev,
  const Array<MultiFab, AMREX_SPACEDIM>& fluxes,
  const MultiFab& U0, MultiFab& Uout)
{
    const auto dx = geom[lev].CellSizeArray();
    const Real dtdx = dt_lev / dx[0];
    const Real dtdy = dt_lev / dx[1];
    for (MFIter mfi(Uout, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u0 = U0.const_array(mfi);
        auto un = Uout.array(mfi);
        auto fx = fluxes[0].const_array(mfi);
        auto fy = fluxes[1].const_array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < det::NCOMP; ++n) {
                  un(i, j, k, n) =
                    u0(i, j, k, n) -
                    dtdx * (fx(i + 1, j, k, n) - fx(i, j, k, n)) -
                    dtdy * (fy(i, j + 1, k, n) - fy(i, j, k, n));
              }
          });
    }
}

void Detonation::RK2Finalize(
  int lev, Real dt_lev,
  const Array<MultiFab, AMREX_SPACEDIM>& fluxes,
  const MultiFab& U0, const MultiFab& Ustage, MultiFab& Uout)
{
    const auto dx = geom[lev].CellSizeArray();
    const Real dtdx = dt_lev / dx[0];
    const Real dtdy = dt_lev / dx[1];
    for (MFIter mfi(Uout, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u0 = U0.const_array(mfi);
        auto us = Ustage.const_array(mfi);
        auto un = Uout.array(mfi);
        auto fx = fluxes[0].const_array(mfi);
        auto fy = fluxes[1].const_array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < det::NCOMP; ++n) {
                  const Real divf =
                    dtdx * (fx(i + 1, j, k, n) - fx(i, j, k, n)) +
                    dtdy * (fy(i, j + 1, k, n) - fy(i, j, k, n));
                  un(i, j, k, n) =
                    0.5_rt * (u0(i, j, k, n) + us(i, j, k, n) - divf);
              }
          });
    }
}

void Detonation::AverageDownFaces(
  Vector<Array<MultiFab, AMREX_SPACEDIM>>& fluxes)
{
    average_faces(fluxes, geom, refRatio(), finest_level);
}

void Detonation::AccumulateFluxRegister(
  const Vector<Array<MultiFab, AMREX_SPACEDIM>>& fluxes,
  Real flux_weight)
{
    if (!do_reflux) {
        return;
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = geom[lev].CellSizeArray();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const Real area = (dir == 0) ? dx[1] : dx[0];
            const Real scale = flux_weight * area;
            if (lev < finest_level && flux_register[lev + 1]) {
                flux_register[lev + 1]->CrseInit(
                  fluxes[lev][dir], dir, 0, 0, det::NCOMP,
                  -scale, FluxRegister::ADD);
            }
            if (lev > 0 && flux_register[lev]) {
                flux_register[lev]->FineAdd(
                  fluxes[lev][dir], dir, 0, 0, det::NCOMP, scale);
            }
        }
    }
}

void Detonation::RefluxAllLevels()
{
    if (!do_reflux) {
        return;
    }
    for (int lev = 0; lev < finest_level; ++lev) {
        if (flux_register[lev + 1]) {
            flux_register[lev + 1]->Reflux(
              state[lev], 1.0_rt, 0, 0, det::NCOMP, geom[lev]);
            flux_register[lev + 1]->setVal(0.0_rt);
        }
    }
}

bool Detonation::CheckFiniteState(const char* stage, Real time)
{
    bool bad = false;
    for (int lev = 0; lev <= finest_level; ++lev) {
        const bool this_level =
          state[lev].contains_nan(0, det::NCOMP, 0, false) ||
          state[lev].contains_inf(0, det::NCOMP, 0, false);
        bad = this_level || bad;
    }
    if (bad) {
        Print() << "NONFINITE_DETECTED stage=" << stage
                << " step=" << istep[0]
                << " t_us=" << time * 1.0e6
                << " finest_level=" << finest_level << '\n';
    }
    return !bad;
}

void Detonation::HydroRK2(Real time, Real dt_lev)
{
    Vector<MultiFab> Ustage(finest_level + 1);
    Vector<MultiFab> Unew(finest_level + 1);
    Vector<Array<MultiFab, AMREX_SPACEDIM>> fluxes(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        Ustage[lev].define(grids[lev], dmap[lev], det::NCOMP, 3);
        Unew[lev].define(grids[lev], dmap[lev], det::NCOMP, 3);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            BoxArray fba = grids[lev];
            fba.surroundingNodes(dir);
            fluxes[lev][dir].define(fba, dmap[lev], det::NCOMP, 0);
        }
    }

    if (do_reflux) {
        for (int lev = 1; lev <= finest_level; ++lev) {
            if (flux_register[lev]) {
                flux_register[lev]->setVal(0.0_rt);
            }
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        ComputeFluxes(lev, time, state, fluxes[lev]);
    }
    AccumulateFluxRegister(fluxes, 0.5_rt * dt_lev);
    for (int lev = 0; lev <= finest_level; ++lev) {
        UpdateWithFluxes(lev, dt_lev, fluxes[lev], state[lev], Ustage[lev]);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        ComputeFluxes(lev, time, Ustage, fluxes[lev]);
    }
    AccumulateFluxRegister(fluxes, 0.5_rt * dt_lev);
    for (int lev = 0; lev <= finest_level; ++lev) {
        RK2Finalize(
          lev, dt_lev, fluxes[lev], state[lev], Ustage[lev], Unew[lev]);
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        state[lev] = std::move(Unew[lev]);
    }
    RefluxAllLevels();
}

void Detonation::ErrorEst(
  int lev, TagBoxArray& tags, Real time, int)
{
    MultiFab sb(grids[lev], dmap[lev], det::NCOMP, 2);
    FillCurrent(lev, time, state, sb, 2);
    const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;
    const Real grad_p_loc = grad_p;
    const Real grad_T_loc = grad_T;
    const Real grad_Y_loc = grad_Y;
    const Real chem_tol_loc = chem_activity_tol;
    const Real front_pre_loc = front_pre;
    const Real front_post_loc = front_post;
    const Real front_loc = front_x_cached;
    const Real zone_lo = refine_zone_lo;
    const Real zone_hi = refine_zone_hi;
    const bool use_zone = (time <= refine_zone_until);
    const Real tmin = Tmin_chem;
    const Real problo_x = Geom(lev).ProbLo(0);
    const Real dx = Geom(lev).CellSize(0);
    Real activity_max = 0.0_rt;
    for (int active_lev = 0; active_lev <= finest_level; ++active_lev) {
        activity_max =
          std::max(activity_max, chem_activity[active_lev].norm0(0, 0, true));
    }
    for (MFIter mfi(sb, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u = sb.const_array(mfi);
        auto chem_act = chem_activity[lev].const_array(mfi);
        auto tag = tags.array(mfi);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real Uloc[det::NCOMP];
              Real Y[detmech::NUM_SPECIES];
              auto primitive = [&](int ii, int jj, Real& pp, Real& TT) {
                  for (int n = 0; n < det::NCOMP; ++n) {
                      Uloc[n] = u(ii, jj, k, n);
                  }
                  Real Tguess = 1000.0_rt;
                  Real rho, vx, vy, cs, gamma;
                  primitive_from_conserved(
                    Uloc, Tguess, rho, vx, vy, pp, TT, cs, gamma, Y);
              };

              Real p0 = 0.0_rt, T0 = 0.0_rt;
              Real pL = 0.0_rt, TL = 0.0_rt;
              Real pR = 0.0_rt, TR = 0.0_rt;
              Real pB = 0.0_rt, TB = 0.0_rt;
              Real pT = 0.0_rt, TT = 0.0_rt;
              primitive(i, j, p0, T0);
              primitive(i - 1, j, pL, TL);
              primitive(i + 1, j, pR, TR);
              primitive(i, j - 1, pB, TB);
              primitive(i, j + 1, pT, TT);

              const Real gpx =
                std::abs(pR - pL) /
                std::max(std::abs(pR) + std::abs(pL), 1.0e-12_rt);
              const Real gpy =
                std::abs(pT - pB) /
                std::max(std::abs(pT) + std::abs(pB), 1.0e-12_rt);
              const Real gtx =
                std::abs(TR - TL) /
                std::max(std::abs(TR) + std::abs(TL), 1.0_rt);
              const Real gty =
                std::abs(TT - TB) /
                std::max(std::abs(TT) + std::abs(TB), 1.0_rt);

              const Real rho =
                std::max(u(i, j, k, det::IRHO), 1.0e-14_rt);
              const Real rhoL =
                std::max(u(i - 1, j, k, det::IRHO), 1.0e-14_rt);
              const Real rhoR =
                std::max(u(i + 1, j, k, det::IRHO), 1.0e-14_rt);
              const Real rhoB =
                std::max(u(i, j - 1, k, det::IRHO), 1.0e-14_rt);
              const Real rhoT =
                std::max(u(i, j + 1, k, det::IRHO), 1.0e-14_rt);
              const Real yh2 = u(i, j, k, det::IFY) / rho;
              const Real yh2L = u(i - 1, j, k, det::IFY) / rhoL;
              const Real yh2R = u(i + 1, j, k, det::IFY) / rhoR;
              const Real yh2B = u(i, j - 1, k, det::IFY) / rhoB;
              const Real yh2T = u(i, j + 1, k, det::IFY) / rhoT;
              const Real yh2o = u(i, j, k, det::IFY + 5) / rho;
              const Real yh2oL = u(i - 1, j, k, det::IFY + 5) / rhoL;
              const Real yh2oR = u(i + 1, j, k, det::IFY + 5) / rhoR;
              const Real yh2oB = u(i, j - 1, k, det::IFY + 5) / rhoB;
              const Real yh2oT = u(i, j + 1, k, det::IFY + 5) / rhoT;
              const Real gyh2x =
                std::abs(yh2R - yh2L) /
                std::max(std::abs(yh2R) + std::abs(yh2L), 1.0e-6_rt);
              const Real gyh2y =
                std::abs(yh2T - yh2B) /
                std::max(std::abs(yh2T) + std::abs(yh2B), 1.0e-6_rt);
              const Real gyh2ox =
                std::abs(yh2oR - yh2oL) /
                std::max(std::abs(yh2oR) + std::abs(yh2oL), 1.0e-6_rt);
              const Real gyh2oy =
                std::abs(yh2oT - yh2oB) /
                std::max(std::abs(yh2oT) + std::abs(yh2oB), 1.0e-6_rt);
              const Real x = problo_x + (i + 0.5_rt) * dx;
              const bool active_reaction =
                (chem_tol_loc > 0.0_rt && activity_max > 0.0_rt &&
                 chem_act(i, j, k) >
                   chem_tol_loc * activity_max);
              const bool front_window =
                (x >= front_loc - front_post_loc &&
                 x <= front_loc + front_pre_loc);
              bool refine =
                (gpx > grad_p_loc || gpy > grad_p_loc ||
                 gtx > grad_T_loc || gty > grad_T_loc ||
                 gyh2x > grad_Y_loc || gyh2y > grad_Y_loc ||
                 gyh2ox > grad_Y_loc || gyh2oy > grad_Y_loc ||
                 active_reaction || front_window);
              if (use_zone) {
                  refine = refine ||
                    (T0 > tmin && x >= zone_lo && x <= zone_hi);
              }
              if (refine) {
                  tag(i, j, k) = tagval;
              } else {
                  tag(i, j, k) = clearval;
              }
          });
    }
    Gpu::streamSynchronize();
}

std::string Detonation::PlotFileName() const
{
    std::ostringstream name;
    name << plot_file << std::fixed << std::setprecision(12) << t_new[0];
    return name.str();
}

void Detonation::WriteStatusFile(
  Real time, Real current_dt, Real wall, Real eta, bool finished)
{
    MultiFab pressure(grids[0], dmap[0], 1, 0);
    MultiFab reaction(grids[0], dmap[0], 1, 0);
    reaction.setVal(0.0_rt);
    for (MFIter mfi(pressure, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto p = pressure.array(mfi);
        auto react = reaction.array(mfi);
        auto u = state[0].const_array(mfi);
        const Real x_ign = x_igniter;
        const Real y_h2_initial = std::max(Y_H2, 1.0e-30_rt);
        const Real problo_x = Geom(0).ProbLo(0);
        const Real dx = Geom(0).CellSize(0);
        const Box reaction_box = amrex::grow(bx, -1);
        ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              Real Uloc[det::NCOMP];
              Real Y[detmech::NUM_SPECIES];
              for (int n = 0; n < det::NCOMP; ++n) {
                  Uloc[n] = u(i, j, k, n);
              }
              Real Tguess = 1000.0_rt;
              Real rho, vx, vy, pp, T, cs, gamma;
              primitive_from_conserved(
                Uloc, Tguess, rho, vx, vy, pp, T, cs, gamma, Y);
              const Real x = problo_x + (i + 0.5_rt) * dx;
              p(i, j, k) = (x >= x_ign) ? pp : 0.0_rt;
          });
        ParallelFor(
          reaction_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              auto progress = [&](int ii, int jj) {
                  const Real rloc =
                    std::max(u(ii, jj, k, det::IRHO), 1.0e-14_rt);
                  const Real yh2 =
                    std::max(u(ii, jj, k, det::IFY), 0.0_rt) / rloc;
                  return std::min(
                    1.0_rt, std::max(0.0_rt, 1.0_rt - yh2 / y_h2_initial));
              };
              const Real gx =
                std::abs(progress(i + 1, j) - progress(i - 1, j));
              const Real gy =
                std::abs(progress(i, j + 1) - progress(i, j - 1));
              const Real x = problo_x + (i + 0.5_rt) * dx;
              react(i, j, k) =
                (x >= x_ign) ? std::max(gx, gy) : 0.0_rt;
          });
        if (reaction_box.isEmpty()) {
            react(bx.smallEnd()) = 0.0_rt;
        }
    }
    const IntVect front_idx = pressure.maxIndex(0, 0);
    const IntVect reaction_idx = reaction.maxIndex(0, 0);
    const Real new_front_x =
      Geom(0).ProbLo(0) + (front_idx[0] + 0.5_rt) * Geom(0).CellSize(0);
    reaction_front_x =
      Geom(0).ProbLo(0) + (reaction_idx[0] + 0.5_rt) * Geom(0).CellSize(0);
    front_separation = std::abs(new_front_x - reaction_front_x);
    if (!finished && time > front_time_previous + 1.0e-30_rt) {
        front_speed =
          (new_front_x - front_x_previous) /
          (time - front_time_previous);
        front_x_previous = new_front_x;
        front_time_previous = time;
    }
    front_x_cached = new_front_x;

    const Real fraction =
      (stop_time > 0.0_rt) ? std::min(1.0_rt, time / stop_time)
                           : (finished ? 1.0_rt : 0.0_rt);
    const Real mass = state[0].sum(det::IRHO);
    const Real energy = state[0].sum(det::IENE);
    Vector<long> level_cells(finest_level + 1);
    long total_cells = 0;
    for (int lev = 0; lev <= finest_level; ++lev) {
        level_cells[lev] = grids[lev].numPts();
        total_cells += level_cells[lev];
    }

    if (ParallelDescriptor::IOProcessor()) {
        const std::string tmp_json = status_file + ".tmp";
        std::ofstream js(tmp_json, std::ios::trunc);
        js << std::setprecision(16)
           << "{\n"
           << "  \"step\": " << istep[0] << ",\n"
           << "  \"t_s\": " << time << ",\n"
           << "  \"t_us\": " << time * 1.0e6 << ",\n"
           << "  \"target_us\": " << stop_time * 1.0e6 << ",\n"
           << "  \"percent\": " << fraction * 100.0_rt << ",\n"
           << "  \"dt_s\": " << current_dt << ",\n"
           << "  \"finest_level\": " << finest_level << ",\n"
           << "  \"cells\": " << total_cells << ",\n"
           << "  \"level_cells\": [";
        for (int lev = 0; lev <= finest_level; ++lev) {
            js << (lev ? ", " : "") << level_cells[lev];
        }
        js << "],\n"
           << "  \"front_x_m\": " << front_x_cached << ",\n"
           << "  \"front_speed_m_s\": " << front_speed << ",\n"
           << "  \"reaction_front_x_m\": " << reaction_front_x << ",\n"
           << "  \"front_separation_mm\": "
           << front_separation * 1.0e3_rt << ",\n"
           << "  \"mass\": " << mass << ",\n"
           << "  \"energy\": " << energy << ",\n"
           << "  \"wall_s\": " << wall << ",\n"
           << "  \"eta_s\": " << eta << ",\n"
           << "  \"finished\": " << (finished ? "true" : "false") << "\n"
           << "}\n";
        js.close();

        std::error_code ec;
        std::filesystem::rename(
          tmp_json, status_file, ec);
        if (ec) {
            std::filesystem::remove(status_file, ec);
            ec.clear();
            std::filesystem::rename(tmp_json, status_file, ec);
        }

        std::ofstream ts(status_text_file, std::ios::trunc);
        ts << "step " << istep[0]
           << "  t " << time * 1.0e6 << " / " << stop_time * 1.0e6
           << " us  " << fraction * 100.0_rt << "%"
           << "  dt " << current_dt << " s"
           << "  finest " << finest_level
           << "  cells " << total_cells
           << "  front " << front_x_cached * 100.0_rt << " cm"
           << "  react " << reaction_front_x * 100.0_rt << " cm"
           << "  gap " << front_separation * 1.0e3_rt << " mm"
           << "  D " << front_speed << " m/s"
           << "  wall " << wall << " s"
           << "  ETA " << eta << " s"
           << (finished ? "  FINISHED" : "") << '\n';
    }
}

void Detonation::WritePlotFile()
{
    constexpr int D_PRESSURE = 0;
    constexpr int D_TEMPERATURE = 1;
    constexpr int D_VELOCITY_X = 2;
    constexpr int D_VELOCITY_Y = 3;
    constexpr int D_DENSITY = 4;
    constexpr int D_GRAD_DENSITY = 5;
    constexpr int D_HEAT_RELEASE = 6;
    constexpr int D_PROGRESS = 7;
    constexpr int D_SPECIES = 8;
    constexpr int NDIAG = D_SPECIES + detmech::NUM_SPECIES;

    const int nlevels = finest_level + 1;
    Vector<MultiFab> diagnostic(nlevels);
    Vector<const MultiFab*> mf(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        diagnostic[lev].define(
          grids[lev], dmap[lev], NDIAG, 0);
        MultiFab sborder(
          grids[lev], dmap[lev], det::NCOMP, 1);
        MultiFab primitive(
          grids[lev], dmap[lev], det::NQ, 1);
        FillCurrent(lev, t_new[0], state, sborder, 1);

        for (MFIter mfi(primitive, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Box& qbox = amrex::grow(bx, 1);
            auto q = primitive.array(mfi);
            auto u = sborder.const_array(mfi);
            ParallelFor(
              qbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  Real Uloc[det::NCOMP];
                  Real Y[detmech::NUM_SPECIES];
                  for (int n = 0; n < det::NCOMP; ++n) {
                      Uloc[n] = u(i, j, k, n);
                  }
                  Real Tguess = 1000.0_rt;
                  Real rho, vx, vy, p, T, cs, gamma;
                  primitive_from_conserved(
                    Uloc, Tguess, rho, vx, vy, p, T, cs, gamma, Y);
                  q(i, j, k, det::QRHO) = rho;
                  q(i, j, k, det::QU) = vx;
                  q(i, j, k, det::QV) = vy;
                  q(i, j, k, det::QPRES) = p;
                  q(i, j, k, det::QTEMP) = T;
                  q(i, j, k, det::QCS) = cs;
                  q(i, j, k, det::QGAM) = gamma;
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      q(i, j, k, det::QY + s) = Y[s];
                  }
              });
        }

        const auto dx = geom[lev].CellSizeArray();
        const Real y_h2_initial = std::max(Y_H2, 1.0e-30_rt);
        for (MFIter mfi(diagnostic[lev], TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto d = diagnostic[lev].array(mfi);
            auto q = primitive.const_array(mfi);
            ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  const Real rho =
                    std::max(q(i, j, k, det::QRHO), 1.0e-14_rt);
                  const Real T = q(i, j, k, det::QTEMP);
                  Real Y[detmech::NUM_SPECIES];
                  Real C[detmech::NUM_SPECIES];
                  Real omega[detmech::NUM_SPECIES];
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      Y[s] = std::max(q(i, j, k, det::QY + s), 0.0_rt);
                  }
                  reaction_rates(rho, Y, T, C, omega);
                  Real heat_release = 0.0_rt;
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      Real es = 0.0_rt;
                      Real cp_s = 0.0_rt;
                      species_thermo(s, T, es, cp_s);
                      const Real hs = es + det::RU * T / detmech::W[s];
                      heat_release -= omega[s] * hs;
                  }
                  const Real drho_dx =
                    0.5_rt *
                    (q(i + 1, j, k, det::QRHO) -
                     q(i - 1, j, k, det::QRHO)) / dx[0];
                  const Real drho_dy =
                    0.5_rt *
                    (q(i, j + 1, k, det::QRHO) -
                     q(i, j - 1, k, det::QRHO)) / dx[1];

                  d(i, j, k, D_PRESSURE) = q(i, j, k, det::QPRES);
                  d(i, j, k, D_TEMPERATURE) = T;
                  d(i, j, k, D_VELOCITY_X) = q(i, j, k, det::QU);
                  d(i, j, k, D_VELOCITY_Y) = q(i, j, k, det::QV);
                  d(i, j, k, D_DENSITY) = rho;
                  d(i, j, k, D_GRAD_DENSITY) =
                    std::sqrt(drho_dx * drho_dx + drho_dy * drho_dy);
                  d(i, j, k, D_HEAT_RELEASE) = heat_release;
                  d(i, j, k, D_PROGRESS) =
                    std::min(1.0_rt, std::max(
                      0.0_rt, 1.0_rt - Y[0] / y_h2_initial));
                  for (int s = 0; s < detmech::NUM_SPECIES; ++s) {
                      d(i, j, k, D_SPECIES + s) = Y[s];
                  }
              });
        }
        mf[lev] = &diagnostic[lev];
    }

    Vector<std::string> names{
      "p", "T", "U_x", "U_y", "rho", "magGradrho", "Qdot", "progress"};
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        names.emplace_back(std::string("Y_") + detmech::species_names[k]);
    }

    const std::filesystem::path plot_path(plot_file);
    const std::filesystem::path output_dir =
      plot_path.has_parent_path() ? plot_path.parent_path()
                                  : std::filesystem::path(".");
    std::error_code ec;
    std::filesystem::create_directories(output_dir, ec);
    WriteMultiLevelPlotfile(
      PlotFileName(), finest_level + 1, mf, names, Geom(),
      t_new[0], istep, refRatio());
    WriteExternalFormats(diagnostic, names);

    if (ParallelDescriptor::IOProcessor()) {
        const std::string catalog = (output_dir / "catalog.json").string();
        if (!std::filesystem::exists(catalog)) {
            std::ofstream js(catalog, std::ios::trunc);
            js << "{\n"
               << "  \"format\": \"AMReX MultiFab plotfile\",\n"
               << "  \"time_unit\": \"s\",\n"
               << "  \"variables\": {\n"
               << "    \"p\": \"Pa\",\n"
               << "    \"T\": \"K\",\n"
               << "    \"U_x\": \"m/s\",\n"
               << "    \"U_y\": \"m/s\",\n"
               << "    \"rho\": \"kg/m^3\",\n"
               << "    \"magGradrho\": \"kg/m^4\",\n"
               << "    \"Qdot\": \"W/m^3\",\n"
               << "    \"progress\": \"-\",\n"
               << "    \"Y_*\": \"mass fraction\"\n"
               << "  },\n"
               << "  \"configuration\": {\n"
               << "    \"stop_time_s\": " << stop_time << ",\n"
               << "    \"max_level\": " << max_level << ",\n"
               << "    \"base_cells_x\": " << geom[0].Domain().length(0)
               << ",\n"
               << "    \"base_cells_y\": " << geom[0].Domain().length(1)
               << ",\n"
               << "    \"cfl\": " << cfl << ",\n"
               << "    \"chem_stride\": " << chem_stride << ",\n"
               << "    \"chem_activity\": " << chem_activity_tol << ",\n"
               << "    \"front_pre_m\": " << front_pre << ",\n"
               << "    \"front_post_m\": " << front_post << ",\n"
               << "    \"export_vtk\": "
               << (export_vtk ? "true" : "false") << ",\n"
               << "    \"export_tecplot\": "
               << (export_tecplot ? "true" : "false") << "\n"
               << "  },\n"
               << "  \"external_formats\": {\n"
               << "    \"paraview_single\": \"paraview/t*.vtm\",\n"
               << "    \"paraview_series\": \"paraview/time_series.pvd\",\n"
               << "    \"paraview_file_series\": \"paraview/fields.series\",\n"
               << "    \"tecplot\": \"tecplot/t*.dat\"\n"
               << "  }\n"
               << "}\n";
        }

        std::string level_text;
        long total_cells = 0;
        for (int lev = 0; lev < nlevels; ++lev) {
            const long cells = grids[lev].numPts();
            total_cells += cells;
            if (lev) {
                level_text += ";";
            }
            level_text += std::to_string(cells);
        }
        const std::string index = (output_dir / "output_index.csv").string();
        const bool write_header = !std::filesystem::exists(index);
        std::ofstream csv(index, std::ios::app);
        std::ostringstream tag;
        tag << 't' << std::fixed << std::setprecision(12) << t_new[0];
        const std::string vtm_file =
          export_vtk
            ? (std::filesystem::path("paraview") /
               (tag.str() + ".vtm")).generic_string()
            : "";
        const std::string tec_file =
          export_tecplot
            ? (std::filesystem::path("tecplot") /
               (tag.str() + ".dat")).generic_string()
            : "";
        if (write_header) {
            csv << "time_s,time_us,step,plotfile,total_cells,level_cells,"
                   "front_x_m,reaction_front_x_m,front_separation_mm,"
                   "mass,energy,paraview_vtm,tecplot_dat\n";
        }
        csv << std::setprecision(16)
            << t_new[0] << ','
            << t_new[0] * 1.0e6 << ','
            << istep[0] << ','
            << std::quoted(PlotFileName()) << ','
            << total_cells << ','
            << std::quoted(level_text) << ','
            << front_x_cached << ','
            << reaction_front_x << ','
            << front_separation * 1.0e3_rt << ','
            << state[0].sum(det::IRHO) << ','
            << state[0].sum(det::IENE) << ','
            << std::quoted(vtm_file) << ','
            << std::quoted(tec_file) << '\n';
    }
}

void Detonation::WriteExternalFormats(
  const Vector<MultiFab>& fields, const Vector<std::string>& names)
{
    const std::filesystem::path plot_path(plot_file);
    const std::filesystem::path output_dir =
      plot_path.has_parent_path() ? plot_path.parent_path()
                                  : std::filesystem::path(".");
    detexport::WriteSnapshot(
      fields, names, Geom(), finest_level, t_new[0],
      output_dir.string(), export_precision,
      export_vtk, export_tecplot);
}

void Detonation::WriteCheckpointFile()
{
    const std::string name = Concatenate(chk_file, istep[0]);
    const int nlevels = finest_level + 1;
    PreBuildDirectorHierarchy(name, "Level_", nlevels, true);
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream hf(name + "/Header");
        hf.precision(17);
        hf << "Detonation checkpoint\n";
        hf << "checkpoint_version " << checkpoint_version << '\n';
        hf << finest_level << '\n';
        for (int s : istep) hf << s << ' ';
        hf << '\n';
        for (Real d : dt) hf << d << ' ';
        hf << '\n';
        for (Real t : t_new) hf << t << ' ';
        hf << '\n';
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(hf);
            hf << '\n';
        }
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Write(
          state[lev], MultiFabFileFullPrefix(lev, name, "Level_", "state"));
    }
}

namespace {
void skip_line(std::istream& is)
{
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}
}

void Detonation::ReadCheckpointFile()
{
    const std::string header = restart_chkfile + "/Header";
    Vector<char> data;
    ParallelDescriptor::ReadAndBcastFile(header, data);
    std::istringstream is(std::string(data.data()));
    std::string line;
    std::getline(is, line);
    std::getline(is, line);
    std::istringstream version_line(line);
    std::string version_tag;
    version_line >> version_tag;
    if (version_tag == "checkpoint_version") {
        int version = 0;
        version_line >> version;
        if (version != checkpoint_version && !allow_legacy_restart) {
            Abort(
              "Checkpoint version " + std::to_string(version) +
              " does not match solver version " +
              std::to_string(checkpoint_version) +
              ". Rerun this case from t=0 with the current executable, or "
              "set det.allow_legacy_restart=1 explicitly.");
        }
        is >> finest_level;
        skip_line(is);
    } else {
        if (!allow_legacy_restart) {
            Abort(
              "Legacy checkpoint has no solver-version metadata. Rerun this "
              "case from t=0 with the current executable, or set "
              "det.allow_legacy_restart=1 explicitly.");
        }
        std::istringstream legacy_line(line);
        legacy_line >> finest_level;
    }
    std::getline(is, line);
    {
        std::istringstream ls(line);
        for (int lev = 0; lev <= finest_level; ++lev) ls >> istep[lev];
    }
    std::getline(is, line);
    {
        std::istringstream ls(line);
        for (int lev = 0; lev <= finest_level; ++lev) ls >> dt[lev];
    }
    std::getline(is, line);
    {
        std::istringstream ls(line);
        for (int lev = 0; lev <= finest_level; ++lev) ls >> t_new[lev];
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray ba;
        ba.readFrom(is);
        skip_line(is);
        DistributionMapping dm(ba, ParallelDescriptor::NProcs());
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);
        state[lev].define(ba, dm, det::NCOMP, 3);
        if (do_reflux && lev > 0) {
            flux_register[lev] = std::make_unique<FluxRegister>(
              ba, dm, refRatio(lev - 1), lev, det::NCOMP);
        }
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(
          state[lev],
          MultiFabFileFullPrefix(
            lev, restart_chkfile, "Level_", "state"));
    }
    for (int lev = 0; lev <= finest_level; ++lev) {
        chem_temperature[lev].define(
          state[lev].boxArray(), state[lev].DistributionMap(), 1, 0);
        chem_hnext[lev].define(
          state[lev].boxArray(), state[lev].DistributionMap(), 1, 0);
        chem_activity[lev].define(
          state[lev].boxArray(), state[lev].DistributionMap(), 1, 0);
        InitializeTemperatureGuess(lev);
    }
}
