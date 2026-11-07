#pragma once

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "MechanismData.H"

#include <cmath>
#include <limits>

namespace det {

using amrex::Real;
using namespace amrex::literals;

static constexpr int IRHO = 0;
static constexpr int IMX = 1;
static constexpr int IMY = 2;
static constexpr int IENE = 3;
static constexpr int IFY = 4;
static constexpr int NCOMP = IFY + detmech::NUM_SPECIES;

static constexpr int QRHO = 0;
static constexpr int QU = 1;
static constexpr int QV = 2;
static constexpr int QPRES = 3;
static constexpr int QTEMP = 4;
static constexpr int QCS = 5;
static constexpr int QGAM = 6;
static constexpr int QY = 7;
static constexpr int NQ = QY + detmech::NUM_SPECIES;
static constexpr int NS = detmech::NUM_SPECIES;
static constexpr int NCHEM = NS + 1;

static constexpr Real RU = Real(8.31446261815324);
static constexpr Real P0 = Real(101325.0);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
species_thermo(int k, Real T, Real& e, Real& cp)
{
    const Real* a = (T <= detmech::Tmid[k]) ? detmech::A_low[k]
                                            : detmech::A_high[k];
    const Real cp_R =
      a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4])));
    const Real h_RT =
      a[0] + T * (a[1] / 2.0_rt + T * (
        a[2] / 3.0_rt + T * (
          a[3] / 4.0_rt + T * a[4] / 5.0_rt))) + a[5] / T;
    cp = RU * cp_R / detmech::W[k];
    e = RU * T * (h_RT - 1.0_rt) / detmech::W[k];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
mixture_thermo(
  const Real* Y, Real T, Real& e, Real& cp, Real& Rmix)
{
    e = 0.0_rt;
    cp = 0.0_rt;
    Rmix = 0.0_rt;
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        Real ek = 0.0_rt;
        Real cpk = 0.0_rt;
        species_thermo(k, T, ek, cpk);
        e += Y[k] * ek;
        cp += Y[k] * cpk;
        Rmix += Y[k] * RU / detmech::W[k];
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real
temperature_from_energy(
  Real e_target, const Real* Y, Real T_guess,
  Real Tmin = 250.0_rt, Real Tmax = 6500.0_rt)
{
    if (!std::isfinite(e_target)) {
        return 0.5_rt * (Tmin + Tmax);
    }

    Real lo = Tmin;
    Real hi = Tmax;
    Real T = std::min(hi, std::max(lo, T_guess));
    for (int it = 0; it < 80; ++it) {
        Real e = 0.0_rt;
        Real cp = 0.0_rt;
        Real Rmix = 0.0_rt;
        mixture_thermo(Y, T, e, cp, Rmix);
        const Real cv = cp - Rmix;
        if (!std::isfinite(e) || !std::isfinite(cv) ||
            cv <= 1.0e-30_rt) {
            return 0.5_rt * (lo + hi);
        }
        const Real residual = e - e_target;
        if (std::abs(residual) <=
            1.0e-11_rt * std::max(std::abs(e_target), 1.0_rt)) {
            break;
        }
        if (residual > 0.0_rt) {
            hi = T;
        } else {
            lo = T;
        }
        Real Tnew = T - residual / cv;
        if (!std::isfinite(Tnew) || Tnew <= lo || Tnew >= hi) {
            Tnew = 0.5_rt * (lo + hi);
        }
        if (std::abs(Tnew - T) <= 1.0e-11_rt * T) {
            T = Tnew;
            break;
        }
        T = Tnew;
    }
    return T;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
primitive_from_conserved(
  const Real* U, Real& Tguess, Real& rho, Real& u, Real& v,
  Real& p, Real& T, Real& cs, Real& gamma, Real* Y)
{
    rho = std::max(U[IRHO], 1.0e-14_rt);
    const Real rhoinv = 1.0_rt / rho;
    u = U[IMX] * rhoinv;
    v = U[IMY] * rhoinv;
    Real sumY = 0.0_rt;
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        Y[k] = std::max(U[IFY + k], 0.0_rt) * rhoinv;
        sumY += Y[k];
    }
    if (sumY <= 1.0e-30_rt) {
        Y[0] = 1.0_rt;
        sumY = 1.0_rt;
    }
    const Real Yinv = 1.0_rt / sumY;
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        Y[k] *= Yinv;
    }

    const Real kinetic = 0.5_rt * rho * (u * u + v * v);
    const Real e = (U[IENE] - kinetic) * rhoinv;
    T = temperature_from_energy(e, Y, Tguess);
    Tguess = T;

    Real cp = 0.0_rt;
    Real Rmix = 0.0_rt;
    Real etmp = 0.0_rt;
    mixture_thermo(Y, T, etmp, cp, Rmix);
    const Real cv = std::max(cp - Rmix, 1.0e-12_rt);
    p = std::max(rho * Rmix * T, 1.0e-12_rt);
    gamma = cp / cv;
    cs = std::sqrt(std::max(gamma * p * rhoinv, 1.0e-12_rt));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real
arrhenius_rate(int r, Real T)
{
    Real k = 0.0_rt;
    for (int t = 0; t < detmech::NTERM; ++t) {
        const Real A = detmech::arr_a[r][t];
        if (A > 0.0_rt) {
            k += A * std::pow(T, detmech::arr_b[r][t]) *
                 std::exp(-detmech::arr_e[r][t] / (RU * T));
        }
    }
    return k;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real
falloff_rate(int r, Real T, Real Mc)
{
    const Real k0 = detmech::k0_a[r] * std::pow(T, detmech::k0_b[r]) *
                    std::exp(-detmech::k0_e[r] / (RU * T));
    const Real ki = detmech::kinf_a[r] * std::pow(T, detmech::kinf_b[r]) *
                    std::exp(-detmech::kinf_e[r] / (RU * T));
    const Real Pr = std::max(k0 * Mc / ki, 1.0e-300_rt);
    Real F = 1.0_rt;
    if (detmech::fall_kind[r] == 2) {
        const Real Fc =
          (1.0_rt - detmech::troe[r][0]) *
            std::exp(-T / detmech::troe[r][1]) +
          detmech::troe[r][0] * std::exp(-T / detmech::troe[r][2]) +
          std::exp(-detmech::troe[r][3] / T);
        const Real lf = std::log10(std::max(Fc, 1.0e-300_rt));
        const Real c = -0.4_rt - 0.67_rt * lf;
        const Real n = 0.75_rt - 1.27_rt * lf;
        const Real lp = std::log10(Pr);
        const Real tmp = (lp + c) / (n - 0.14_rt * (lp + c));
        F = std::pow(10.0_rt, lf / (1.0_rt + tmp * tmp));
    }
    return ki * (Pr / (1.0_rt + Pr)) * F;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real
falloff_rate_parts(
  int r, Real T, Real Mc, Real& Pr, Real& F, Real& k0)
{
    k0 = detmech::k0_a[r] * std::pow(T, detmech::k0_b[r]) *
         std::exp(-detmech::k0_e[r] / (RU * T));
    const Real ki = detmech::kinf_a[r] * std::pow(T, detmech::kinf_b[r]) *
                    std::exp(-detmech::kinf_e[r] / (RU * T));
    Pr = std::max(k0 * Mc / ki, 1.0e-300_rt);
    F = 1.0_rt;
    if (detmech::fall_kind[r] == 2) {
        const Real Fc =
          (1.0_rt - detmech::troe[r][0]) *
            std::exp(-T / detmech::troe[r][1]) +
          detmech::troe[r][0] * std::exp(-T / detmech::troe[r][2]) +
          std::exp(-detmech::troe[r][3] / T);
        const Real lf = std::log10(std::max(Fc, 1.0e-300_rt));
        const Real c = -0.4_rt - 0.67_rt * lf;
        const Real n = 0.75_rt - 1.27_rt * lf;
        const Real lp = std::log10(Pr);
        const Real tmp = (lp + c) / (n - 0.14_rt * (lp + c));
        F = std::pow(10.0_rt, lf / (1.0_rt + tmp * tmp));
    }
    return ki * (Pr / (1.0_rt + Pr)) * F;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
reaction_rates_jac(
  Real T, const Real* C, Real* omega, Real Jac[NS][NS],
  Real* dTvec, Real* u, Real& cv_vol)
{
    Real h[NS];
    Real s[NS];
    Real cp[NS];
    for (int k = 0; k < NS; ++k) {
        const Real* a = (T <= detmech::Tmid[k]) ? detmech::A_low[k]
                                                : detmech::A_high[k];
        cp[k] = RU * (a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4]))));
        h[k] = RU * (T * (a[0] + T * (a[1] / 2.0_rt + T * (
          a[2] / 3.0_rt + T * (a[3] / 4.0_rt + T * a[4] / 5.0_rt)))) + a[5]);
        s[k] = RU * (a[0] * std::log(T) + T * (
          a[1] + T * (a[2] / 2.0_rt + T * (
            a[3] / 3.0_rt + T * a[4] / 4.0_rt))) + a[6]);
        u[k] = h[k] - RU * T;
        omega[k] = 0.0_rt;
        dTvec[k] = 0.0_rt;
        for (int j = 0; j < NS; ++j) {
            Jac[k][j] = 0.0_rt;
        }
    }
    cv_vol = 0.0_rt;
    for (int k = 0; k < NS; ++k) {
        cv_vol += C[k] * (cp[k] - RU);
    }

    const Real invT = 1.0_rt / T;
    const Real ln_c0 = std::log(P0 / (RU * T));
    Real pder[2 * NS];
    for (int r = 0; r < detmech::NUM_REACTIONS; ++r) {
        Real Mc = 0.0_rt;
        if (detmech::nbody[r] > 0) {
            for (int k = 0; k < NS; ++k) {
                Mc += detmech::eff[r][k] * C[k];
            }
        }

        Real kf = 0.0_rt;
        Real dlnkf = 0.0_rt;
        Real dkfdM = 0.0_rt;
        if (detmech::nbody[r] == 2) {
            Real Pr = 0.0_rt, F = 0.0_rt, k0 = 0.0_rt;
            kf = falloff_rate_parts(r, T, Mc, Pr, F, k0);
            const Real Pr_T = Pr;
            const Real F_T = F;
            const Real k0_T = k0;
            const Real kp =
              falloff_rate_parts(r, 1.0000005_rt * T, Mc, Pr, F, k0);
            const Real km =
              falloff_rate_parts(r, 0.9999995_rt * T, Mc, Pr, F, k0);
            dlnkf = (std::log(std::max(kp, 1.0e-300_rt)) -
                     std::log(std::max(km, 1.0e-300_rt))) /
                    (1.0e-6_rt * T);
            const Real one_plus = 1.0_rt + Pr_T;
            dkfdM = F_T * k0_T / (one_plus * one_plus);
        } else {
            kf = arrhenius_rate(r, T);
            Real dkdt = 0.0_rt;
            for (int t = 0; t < detmech::NTERM; ++t) {
                const Real A = detmech::arr_a[r][t];
                if (A > 0.0_rt) {
                    const Real b = detmech::arr_b[r][t];
                    const Real E = detmech::arr_e[r][t];
                    dkdt += A * std::pow(T, b) *
                            std::exp(-E / (RU * T)) *
                            (b / T + E / (RU * T * T));
                }
            }
            dlnkf = (kf > 0.0_rt) ? dkdt / kf : 0.0_rt;
        }

        Real pf = 1.0_rt;
        for (int a = 0; a < detmech::NR_MAX; ++a) {
            const int k = detmech::r_idx[r][a];
            if (k < 0) {
                break;
            }
            pf *= std::pow(std::max(C[k], 0.0_rt), detmech::r_nu[r][a]);
        }
        Real pr = 1.0_rt;
        for (int a = 0; a < detmech::NP_MAX; ++a) {
            const int k = detmech::p_idx[r][a];
            if (k < 0) {
                break;
            }
            pr *= std::pow(std::max(C[k], 0.0_rt), detmech::p_nu[r][a]);
        }

        Real qf = kf * pf;
        Real qr = 0.0_rt;
        if (detmech::nbody[r] == 1) {
            qf *= Mc;
        }
        Real dH = 0.0_rt;
        Real dS = 0.0_rt;
        Real dn = 0.0_rt;
        for (int a = 0; a < detmech::NN_MAX; ++a) {
            const int k = detmech::n_idx[r][a];
            if (k < 0) {
                break;
            }
            const Real nu = detmech::n_nu[r][a];
            dH += nu * h[k];
            dS += nu * s[k];
            dn += nu;
        }
        Real lnKc = dS / RU - dH * invT / RU + dn * ln_c0;
        lnKc = std::min(150.0_rt, std::max(-150.0_rt, lnKc));
        const Real kr = kf * std::exp(-lnKc) * detmech::revflag[r];
        qr = kr * pr;
        if (detmech::nbody[r] == 1) {
            qr *= Mc;
        }
        const Real dlnKc = dH * invT * invT / RU - dn * invT;
        const Real dlnkr = dlnkf - dlnKc;

        for (int a = 0; a < detmech::NN_MAX; ++a) {
            const int k = detmech::n_idx[r][a];
            if (k < 0) {
                break;
            }
            const Real nu = detmech::n_nu[r][a];
            omega[k] += nu * (qf - qr);
            dTvec[k] += nu * (qf * dlnkf - qr * dlnkr);
        }

        for (int j = 0; j < 2 * NS; ++j) {
            pder[j] = 0.0_rt;
        }
        for (int a = 0; a < detmech::NR_MAX; ++a) {
            const int k = detmech::r_idx[r][a];
            if (k < 0) {
                break;
            }
            if (C[k] > 0.0_rt) {
                pder[k] += qf * detmech::r_nu[r][a] / C[k];
            }
        }
        if (detmech::nbody[r] == 1) {
            for (int k = 0; k < NS; ++k) {
                pder[k] += kf * pf * detmech::eff[r][k];
            }
        } else if (detmech::nbody[r] == 2) {
            for (int k = 0; k < NS; ++k) {
                pder[k] += dkfdM * pf * detmech::eff[r][k];
            }
        }
        for (int a = 0; a < detmech::NP_MAX; ++a) {
            const int k = detmech::p_idx[r][a];
            if (k < 0) {
                break;
            }
            if (C[k] > 0.0_rt) {
                pder[NS + k] += qr * detmech::p_nu[r][a] / C[k];
            }
        }
        if (detmech::nbody[r] == 1) {
            for (int k = 0; k < NS; ++k) {
                pder[NS + k] += kr * pr * detmech::eff[r][k];
            }
        } else if (detmech::nbody[r] == 2) {
            for (int k = 0; k < NS; ++k) {
                pder[NS + k] += dkfdM * pr * detmech::eff[r][k];
            }
        }
        for (int a = 0; a < detmech::NN_MAX; ++a) {
            const int k = detmech::n_idx[r][a];
            if (k < 0) {
                break;
            }
            const Real nu = detmech::n_nu[r][a];
            for (int j = 0; j < NS; ++j) {
                Jac[k][j] += nu * (pder[j] - pder[NS + j]);
            }
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
reaction_rates(
  Real rho, const Real* Y, Real T, Real* C, Real* omega)
{
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        C[k] = rho * Y[k] / detmech::W[k];
        omega[k] = 0.0_rt;
    }

    Real h[detmech::NUM_SPECIES];
    Real s[detmech::NUM_SPECIES];
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        const Real* a = (T <= detmech::Tmid[k]) ? detmech::A_low[k]
                                                : detmech::A_high[k];
        h[k] = RU * (T * (a[0] + T * (a[1] / 2.0_rt + T * (
          a[2] / 3.0_rt + T * (a[3] / 4.0_rt + T * a[4] / 5.0_rt)))) + a[5]);
        s[k] = RU * (a[0] * std::log(T) + T * (
          a[1] + T * (a[2] / 2.0_rt + T * (
            a[3] / 3.0_rt + T * a[4] / 4.0_rt))) + a[6]);
    }

    for (int r = 0; r < detmech::NUM_REACTIONS; ++r) {
        Real Mc = 0.0_rt;
        if (detmech::nbody[r] > 0) {
            for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
                Mc += detmech::eff[r][k] * C[k];
            }
        }
        const Real kf_rate = (detmech::nbody[r] == 2)
                               ? falloff_rate(r, T, Mc)
                               : arrhenius_rate(r, T);
        Real qf = kf_rate;
        if (detmech::nbody[r] == 1) {
            qf *= Mc;
        }
        for (int a = 0; a < detmech::NR_MAX; ++a) {
            const int k = detmech::r_idx[r][a];
            if (k < 0) {
                break;
            }
            qf *= std::pow(std::max(C[k], 0.0_rt),
                           detmech::r_nu[r][a]);
        }

        Real dH = 0.0_rt;
        Real dS = 0.0_rt;
        Real dn = 0.0_rt;
        for (int a = 0; a < detmech::NN_MAX; ++a) {
            const int k = detmech::n_idx[r][a];
            if (k < 0) {
                break;
            }
            const Real nu = detmech::n_nu[r][a];
            dH += nu * h[k];
            dS += nu * s[k];
            dn += nu;
        }
        Real lnKc = dS / RU - dH / (RU * T) +
                    dn * std::log(P0 / (RU * T));
        lnKc = std::min(150.0_rt, std::max(-150.0_rt, lnKc));
        const Real kr_rate =
          kf_rate * std::exp(-lnKc) * detmech::revflag[r];
        Real qr = kr_rate;
        if (detmech::nbody[r] == 1) {
            qr *= Mc;
        }
        for (int a = 0; a < detmech::NP_MAX; ++a) {
            const int k = detmech::p_idx[r][a];
            if (k < 0) {
                break;
            }
            qr *= std::pow(std::max(C[k], 0.0_rt),
                           detmech::p_nu[r][a]);
        }

        for (int a = 0; a < detmech::NN_MAX; ++a) {
            const int k = detmech::n_idx[r][a];
            if (k < 0) {
                break;
            }
            omega[k] += detmech::n_nu[r][a] * (qf - qr);
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool
solve_linear_system(Real A[NCHEM][NCHEM], Real* b, Real* x)
{
    for (int col = 0; col < NCHEM; ++col) {
        int piv = col;
        Real amax = std::abs(A[col][col]);
        for (int rr = col + 1; rr < NCHEM; ++rr) {
            if (std::abs(A[rr][col]) > amax) {
                amax = std::abs(A[rr][col]);
                piv = rr;
            }
        }
        if (amax < 1.0e-300_rt) {
            return false;
        }
        if (piv != col) {
            for (int cc = col; cc < NCHEM; ++cc) {
                const Real tmp = A[col][cc];
                A[col][cc] = A[piv][cc];
                A[piv][cc] = tmp;
            }
            const Real tmp = b[col];
            b[col] = b[piv];
            b[piv] = tmp;
        }
        for (int rr = col + 1; rr < NCHEM; ++rr) {
            const Real fac = A[rr][col] / A[col][col];
            if (fac != 0.0_rt) {
                for (int cc = col; cc < NCHEM; ++cc) {
                    A[rr][cc] -= fac * A[col][cc];
                }
                b[rr] -= fac * b[col];
            }
        }
    }
    for (int rr = NCHEM - 1; rr >= 0; --rr) {
        Real sum = b[rr];
        for (int cc = rr + 1; cc < NCHEM; ++cc) {
            sum -= A[rr][cc] * x[cc];
        }
        if (std::abs(A[rr][rr]) < 1.0e-300_rt) {
            return false;
        }
        x[rr] = sum / A[rr][rr];
    }
    return true;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool
implicit_chemistry_step(
  const Real* C0, Real rhoe0, Real T0, Real dt,
  Real tol_y, Real tol_T, Real Tmin, Real Tmax,
  Real* Cout, Real& Tout, Real& change_y, Real& change_T,
  Real* debug = nullptr)
{
    Real x[NCHEM];
    Real xnew[NCHEM];
    Real G[NCHEM];
    Real rhs[NCHEM];
    Real dx[NCHEM];
    Real J[NCHEM][NCHEM];
    Real om[NS];
    Real dTvec[NS];
    Real u[NS];
    Real jac[NS][NS];

    Real Ctot = 0.0_rt;
    for (int k = 0; k < NS; ++k) {
        x[k] = C0[k];
        Ctot += C0[k];
    }
    x[NS] = std::min(Tmax, std::max(Tmin, T0));
    if (Ctot < 1.0e-12_rt) {
        return false;
    }
    const Real Cinv = 1.0_rt / Ctot;
    const Real Einv = 1.0_rt / std::max(std::abs(rhoe0), 1.0_rt);
    Real trial_Y[NS];
    Real trial_C[NS];
    Real trial_omega[NS];

    auto evaluate_residual =
      [&](const Real* C, Real T) {
          Real rho = 0.0_rt;
          for (int k = 0; k < NS; ++k) {
              rho += C[k] * detmech::W[k];
          }
          if (!(rho > 0.0_rt) || !std::isfinite(rho)) {
              return std::numeric_limits<Real>::max();
          }
          for (int k = 0; k < NS; ++k) {
              trial_Y[k] = C[k] * detmech::W[k] / rho;
          }
          reaction_rates(rho, trial_Y, T, trial_C, trial_omega);
          Real e_vol = 0.0_rt;
          for (int k = 0; k < NS; ++k) {
              Real e_k = 0.0_rt;
              Real cp_k = 0.0_rt;
              species_thermo(k, T, e_k, cp_k);
              e_vol += C[k] * e_k * detmech::W[k];
          }
          Real next_residual = 0.0_rt;
          for (int k = 0; k < NS; ++k) {
              const Real g = (C[k] - C0[k] - dt * trial_omega[k]) * Cinv;
              next_residual += g * g;
          }
          const Real ge = (e_vol - rhoe0) * Einv;
          next_residual += ge * ge;
          return next_residual;
      };

    for (int it = 0; it < 25; ++it) {
        Real cv_vol = 0.0_rt;
        reaction_rates_jac(x[NS], x, om, jac, dTvec, u, cv_vol);

        Real e_vol = 0.0_rt;
        for (int k = 0; k < NS; ++k) {
            G[k] = (x[k] - C0[k] - dt * om[k]) * Cinv;
            e_vol += x[k] * u[k];
            for (int j = 0; j < NS; ++j) {
                J[k][j] =
                  ((k == j ? 1.0_rt : 0.0_rt) - dt * jac[k][j]) * Cinv;
            }
            J[k][NS] = -dt * dTvec[k] * Cinv;
        }
        G[NS] = (e_vol - rhoe0) * Einv;
        for (int j = 0; j < NS; ++j) {
            J[NS][j] = u[j] * Einv;
        }
        J[NS][NS] = cv_vol * Einv;

        Real residual = 0.0_rt;
        for (int k = 0; k < NCHEM; ++k) {
            residual += G[k] * G[k];
            rhs[k] = -G[k];
            dx[k] = 0.0_rt;
        }
        if (!std::isfinite(residual)) {
            if (debug) {
                debug[0] = it;
                debug[1] = residual;
            }
            return false;
        }
        if (debug) {
            debug[0] = it;
            debug[1] = residual;
        }
        if (residual < 1.0e-12_rt) {
            break;
        }
        if (!solve_linear_system(J, rhs, dx)) {
            if (debug) {
                debug[2] = -1.0_rt;
            }
            return false;
        }

        Real lambda = 1.0_rt;
        bool accepted = false;
        for (int ls = 0; ls < 10; ++ls) {
            for (int k = 0; k < NS; ++k) {
                xnew[k] = std::max(0.0_rt, x[k] + lambda * dx[k]);
            }
            xnew[NS] =
              std::min(Tmax, std::max(Tmin, x[NS] + lambda * dx[NS]));

            const Real next_residual = evaluate_residual(xnew, xnew[NS]);
            if (std::isfinite(next_residual) && next_residual < residual) {
                accepted = true;
                break;
            }
            lambda *= 0.5_rt;
        }
        if (!accepted) {
            if (debug) {
                debug[2] = -2.0_rt;
            }
            return false;
        }
        for (int k = 0; k < NCHEM; ++k) {
            x[k] = xnew[k];
        }
    }

    const Real residual = evaluate_residual(x, x[NS]);
    if (!(residual < 1.0e-12_rt)) {
        return false;
    }
    for (int k = 0; k < NS; ++k) {
        if (!std::isfinite(x[k]) || x[k] < 0.0_rt) {
            return false;
        }
        Cout[k] = x[k];
    }
    if (!std::isfinite(x[NS]) || x[NS] < Tmin || x[NS] > Tmax) {
        return false;
    }
    Tout = x[NS];
    change_y = 0.0_rt;
    for (int k = 0; k < NS; ++k) {
        change_y = std::max(change_y, std::abs(Cout[k] - C0[k]) * Cinv);
    }
    change_T = std::abs(Tout - T0) / std::max(T0, 1.0_rt);
    return true;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
chemistry_step(
  Real rho, Real e0, Real dt, Real Tmin, Real hmax,
  Real tol_y, Real tol_T, Real* Y, Real& T, Real& hnext,
  Real& activity)
{
    activity = 0.0_rt;
    if (dt <= 0.0_rt) {
        return;
    }
    if (T < Tmin) {
        hnext = hmax;
        return;
    }

    Real h = std::min(
      hmax, std::max(1.0e-15_rt,
        std::isfinite(hnext) ? hnext : hmax));
    Real time = 0.0_rt;
    int nsub = 0;
    while (time < dt - 1.0e-300_rt && nsub < 64) {
        Real h_use = std::min(h, dt - time);
        bool truncated = (h_use < 0.999_rt * h);
        const Real Tstep = T;
        Real C0[NS];
        Real Cnew[NS];
        for (int k = 0; k < NS; ++k) {
            C0[k] = rho * std::max(Y[k], 0.0_rt) / detmech::W[k];
        }
        bool ok = false;
        Real Tnew = T;
        Real dy = 0.0_rt;
        Real dT = 0.0_rt;
        for (int attempt = 0; attempt < 12; ++attempt) {
            ok = implicit_chemistry_step(
              C0, rho * e0, T, h_use, tol_y, tol_T, 250.0_rt, 6500.0_rt,
              Cnew, Tnew, dy, dT, nullptr);
            if (ok && dy <= tol_y && dT <= tol_T) {
                break;
            }
            h_use *= 0.5_rt;
            h = h_use;
            truncated = false;
            if (h_use < 1.0e-16_rt * std::max(dt, 1.0e-30_rt)) {
                break;
            }
        }
        if (!ok) {
            break;
        }
        {
            Real sum_y = 0.0_rt;
            for (int k = 0; k < NS; ++k) {
                Y[k] = Cnew[k] * detmech::W[k] / rho;
                sum_y += Y[k];
            }
            if (sum_y > 0.0_rt) {
                for (int k = 0; k < NS; ++k) {
                    Y[k] /= sum_y;
                }
            }
            T = Tnew;
        }
        activity =
          std::max(activity, std::abs(T - Tstep) / std::max(Tstep, 1.0_rt));
        time += h_use;
        ++nsub;
        const Real fy =
          std::sqrt(tol_y / std::max(dy, 1.0e-12_rt));
        const Real fT =
          std::sqrt(tol_T / std::max(dT, 1.0e-12_rt));
        Real fac = 0.9_rt * std::min(fy, fT);
        fac = std::min(2.0_rt, std::max(0.5_rt, fac));
        if (!truncated) {
            h = std::min(hmax, h_use * fac);
        }
    }
    hnext = std::min(hmax, std::max(h, 1.0e-15_rt));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
hllc_flux(
  int normal, Real cL, Real cR,
  Real rhoL, Real uL, Real vL, Real pL, Real EL, const Real* YL,
  Real rhoR, Real uR, Real vR, Real pR, Real ER, const Real* YR,
  Real* F)
{
    Real unL = (normal == 0) ? uL : vL;
    Real utL = (normal == 0) ? vL : uL;
    Real unR = (normal == 0) ? uR : vR;
    Real utR = (normal == 0) ? vR : uR;
    const Real sqL = std::sqrt(rhoL);
    const Real sqR = std::sqrt(rhoR);
    const Real ua = (sqL * unL + sqR * unR) / (sqL + sqR);
    const Real ca = (sqL * cL + sqR * cR) / (sqL + sqR);
    const Real SL = std::min(unL - cL, ua - ca);
    const Real SR = std::max(unR + cR, ua + ca);

    Real FL[det::NCOMP];
    Real FR[det::NCOMP];
    if (normal == 0) {
        FL[IRHO] = rhoL * unL;
        FL[IMX] = rhoL * unL * unL + pL;
        FL[IMY] = rhoL * unL * utL;
        FL[IENE] = unL * (EL + pL);
        FR[IRHO] = rhoR * unR;
        FR[IMX] = rhoR * unR * unR + pR;
        FR[IMY] = rhoR * unR * utR;
        FR[IENE] = unR * (ER + pR);
    } else {
        FL[IRHO] = rhoL * unL;
        FL[IMX] = rhoL * unL * utL;
        FL[IMY] = rhoL * unL * unL + pL;
        FL[IENE] = unL * (EL + pL);
        FR[IRHO] = rhoR * unR;
        FR[IMX] = rhoR * unR * utR;
        FR[IMY] = rhoR * unR * unR + pR;
        FR[IENE] = unR * (ER + pR);
    }
    for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
        FL[IFY + k] = rhoL * unL * YL[k];
        FR[IFY + k] = rhoR * unR * YR[k];
    }

    auto set_rusanov = [&]() {
        Real UL[NCOMP];
        Real UR[NCOMP];
        UL[IRHO] = rhoL;
        UR[IRHO] = rhoR;
        if (normal == 0) {
            UL[IMX] = rhoL * unL;
            UL[IMY] = rhoL * utL;
            UR[IMX] = rhoR * unR;
            UR[IMY] = rhoR * utR;
        } else {
            UL[IMX] = rhoL * utL;
            UL[IMY] = rhoL * unL;
            UR[IMX] = rhoR * utR;
            UR[IMY] = rhoR * unR;
        }
        UL[IENE] = EL;
        UR[IENE] = ER;
        for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
            UL[IFY + k] = rhoL * YL[k];
            UR[IFY + k] = rhoR * YR[k];
        }
        const Real amax = std::max(
          std::abs(unL) + cL, std::abs(unR) + cR);
        for (int n = 0; n < NCOMP; ++n) {
            F[n] = 0.5_rt * (FL[n] + FR[n] -
                              amax * (UR[n] - UL[n]));
        }
    };

    auto finish_flux = [&]() {
        bool bad = !std::isfinite(SL) || !std::isfinite(SR);
        for (int n = 0; n < NCOMP; ++n) {
            bad = bad || !std::isfinite(F[n]);
        }
        if (bad) {
            set_rusanov();
        }
    };

    if (SL >= 0.0_rt) {
        for (int n = 0; n < NCOMP; ++n) {
            F[n] = FL[n];
        }
        finish_flux();
        return;
    }
    if (SR <= 0.0_rt) {
        for (int n = 0; n < NCOMP; ++n) {
            F[n] = FR[n];
        }
        finish_flux();
        return;
    }
    const Real den = rhoL * (SL - unL) - rhoR * (SR - unR);
    if (!std::isfinite(den) || std::abs(den) < 1.0e-30_rt) {
        set_rusanov();
        return;
    }
    const Real SM =
      (pR - pL + rhoL * unL * (SL - unL) - rhoR * unR * (SR - unR)) /
      den;
    if (SM >= 0.0_rt) {
        const Real fac = rhoL * (SL - unL) / (SL - SM);
        const Real us = fac * SM;
        const Real ut = fac * utL;
        const Real Us_rho = fac;
        const Real Us_mx = (normal == 0) ? us : ut;
        const Real Us_my = (normal == 0) ? ut : us;
        const Real Us_E = fac * (
          EL / rhoL + (SM - unL) * (SM + pL / (rhoL * (SL - unL))));
        F[IRHO] = FL[IRHO] + SL * (Us_rho - rhoL);
        F[IMX] = FL[IMX] + SL * (Us_mx - rhoL * uL);
        F[IMY] = FL[IMY] + SL * (Us_my - rhoL * vL);
        F[IENE] = FL[IENE] + SL * (Us_E - EL);
        for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
            F[IFY + k] = FL[IFY + k] +
                         SL * (fac * YL[k] - rhoL * YL[k]);
        }
    } else {
        const Real fac = rhoR * (SR - unR) / (SR - SM);
        const Real us = fac * SM;
        const Real ut = fac * utR;
        const Real Us_rho = fac;
        const Real Us_mx = (normal == 0) ? us : ut;
        const Real Us_my = (normal == 0) ? ut : us;
        const Real Us_E = fac * (
          ER / rhoR + (SM - unR) * (SM + pR / (rhoR * (SR - unR))));
        F[IRHO] = FR[IRHO] + SR * (Us_rho - rhoR);
        F[IMX] = FR[IMX] + SR * (Us_mx - rhoR * uR);
        F[IMY] = FR[IMY] + SR * (Us_my - rhoR * vR);
        F[IENE] = FR[IENE] + SR * (Us_E - ER);
        for (int k = 0; k < detmech::NUM_SPECIES; ++k) {
            F[IFY + k] = FR[IFY + k] +
                         SR * (fac * YR[k] - rhoR * YR[k]);
        }
    }
    finish_flux();
}

} // namespace det
