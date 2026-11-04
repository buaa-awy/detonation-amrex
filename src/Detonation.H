#pragma once

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_FillPatcher.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

#include "Physics.H"

class Detonation final : public amrex::AmrCore
{
public:
    Detonation();
    ~Detonation() override = default;

    void Evolve();
    void InitData();

    void MakeNewLevelFromCoarse(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override;
    void RemakeLevel(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override;
    void ClearLevel(int lev) override;
    void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override;
    void ErrorEst(
      int lev, amrex::TagBoxArray& tags, amrex::Real time,
      int ngrow) override;

public:
    void ReadParameters();
    void InitLevelData(int lev);
    void InitializeTemperatureGuess(int lev);
    void UpdateFineMasks();
    void FillCurrent(
      int lev, amrex::Real time,
      const amrex::Vector<amrex::MultiFab>& source,
      amrex::MultiFab& data, int ngrow);
    void FillCoarsePatch(
      int lev, amrex::Real time, amrex::MultiFab& data);

    void AverageDown();
    amrex::Real ComputeDt();
    amrex::Real EstTimeStep(int lev);

    void ChemistryHalfStep(amrex::Real dt);
    void HydroRK2(amrex::Real time, amrex::Real dt);
    void ComputeFluxes(
      int lev, amrex::Real time,
      const amrex::Vector<amrex::MultiFab>& source,
      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fluxes);
    void UpdateWithFluxes(
      int lev, amrex::Real dt,
      const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fluxes,
      const amrex::MultiFab& U0, amrex::MultiFab& Uout);
    void RK2Finalize(
      int lev, amrex::Real dt,
      const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fluxes,
      const amrex::MultiFab& U0, const amrex::MultiFab& Ustage,
      amrex::MultiFab& Uout);
    void AverageDownFaces(
      amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>& fluxes);
    void AccumulateFluxRegister(
      const amrex::Vector<
        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>& fluxes,
      amrex::Real flux_weight);
    void RefluxAllLevels();
    bool CheckFiniteState(const char* stage, amrex::Real time);

    void WritePlotFile();
    void WriteExternalFormats(
      const amrex::Vector<amrex::MultiFab>& fields,
      const amrex::Vector<std::string>& names);
    void WriteCheckpointFile();
    void ReadCheckpointFile();
    void WriteStatusFile(
      amrex::Real time, amrex::Real current_dt, amrex::Real wall,
      amrex::Real eta, bool finished);
    amrex::Real EstimateEta(
      amrex::Real time, amrex::Real wall, bool finished);
    std::string PlotFileName() const;

private:
    amrex::Vector<amrex::MultiFab> state;
    amrex::Vector<amrex::MultiFab> chem_temperature;
    amrex::Vector<amrex::MultiFab> chem_hnext;
    amrex::Vector<amrex::MultiFab> chem_activity;
    amrex::Vector<amrex::iMultiFab> fine_mask;
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;
    amrex::Vector<int> istep;
    amrex::Vector<amrex::BCRec> bcs;
    amrex::Vector<std::unique_ptr<amrex::FillPatcher<amrex::MultiFab>>>
      fillpatcher;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> flux_register;

    int max_step = 1000000;
    int regrid_int = 20;
    int plot_int = 100;
    int chk_int = 100;
    amrex::Real plot_time_int = 0.0;
    amrex::Real chk_time_int = 0.0;
    int progress_int = 1;
    int status_int = 1;
    int chem_stride = 1;
    int export_precision = 10;
    bool export_vtk = true;
    bool export_tecplot = true;
    bool do_reflux = true;
    bool verbose = false;
    amrex::Real stop_time = 1.0e-4;
    amrex::Real next_plot_time = 0.0;
    amrex::Real next_chk_time = 0.0;
    amrex::Real last_plot_time = -1.0;
    amrex::Real last_chk_time = -1.0;
    amrex::Real chem_block_dt = 0.0;
    amrex::Real eta_smoothed = -1.0;
    amrex::Real eta_last_time = -1.0;
    amrex::Real eta_last_wall = -1.0;
    amrex::Real cfl = 0.3;
    amrex::Real change_max = 1.1;
    amrex::Real wall_start = 0.0;
    std::string plot_file = "plt";
    std::string chk_file = "chk";
    std::string status_file = "status.json";
    std::string status_text_file = "status.txt";
    std::string restart_chkfile;

    amrex::Real x_igniter = 0.01;
    amrex::Real T_igniter = 1800.0;
    amrex::Real p_igniter = 101325.0;
    amrex::Real T_fuel = 300.0;
    amrex::Real p_fuel = 8000.0;
    amrex::Real perturb_amp = 0.12;
    int perturb_waves = 4;
    amrex::Real Y_H2 = 0.1119067440;
    amrex::Real Y_O2 = 0.8880932560;
    amrex::Real grad_p = 0.15;
    amrex::Real grad_T = 0.10;
    amrex::Real grad_Y = 0.02;
    amrex::Real refine_zone_lo = 0.009;
    amrex::Real refine_zone_hi = 0.013;
    amrex::Real refine_zone_until = 5.0e-6;
    amrex::Real front_regrid_distance = 5.0e-4;
    amrex::Real last_regrid_front = -1.0;
    amrex::Real front_x_previous = 0.01;
    amrex::Real front_time_previous = 0.0;
    amrex::Real front_speed = 0.0;
    amrex::Real reaction_front_x = 0.0;
    amrex::Real front_separation = 0.0;
    amrex::Real chem_activity_tol = 0.05;
    amrex::Real front_pre = 1.5e-3;
    amrex::Real front_post = 2.5e-3;
    amrex::Real Tmin_chem = 650.0;
    amrex::Real hmax_chem = 1.0e-8;
    amrex::Real tol_y_chem = 0.05;
    amrex::Real tol_T_chem = 0.03;
    amrex::Real front_x_cached = 0.01;
    bool allow_legacy_restart = false;

    static constexpr int checkpoint_version = 2;
};
