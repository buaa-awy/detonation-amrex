#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <iomanip>
#include <iostream>
#include <vector>

#include "Physics.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr int nout = det::NS + 11;
        amrex::Gpu::DeviceVector<amrex::Real> device(nout);
        amrex::Gpu::streamSynchronize();
        auto* out = device.data();
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
            using namespace det;
            Real Y[NS] = {0.0_rt};
            Y[0] = 0.1119067440_rt;
            Y[3] = 0.8880932560_rt;
            const Real T0 = 1800.0_rt;
            const Real p0 = 151987.5_rt;
            Real sum_inv_w = 0.0_rt;
            for (int k = 0; k < NS; ++k) {
                sum_inv_w += Y[k] / detmech::W[k];
            }
            const Real rho = p0 / (RU * T0) / sum_inv_w;
            Real e0 = 0.0_rt;
            Real cp = 0.0_rt;
            Real Rmix = 0.0_rt;
            mixture_thermo(Y, T0, e0, cp, Rmix);

            Real C0[NS];
            Real C1[NS];
            Real om0[NS];
            for (int k = 0; k < NS; ++k) {
                C0[k] = rho * Y[k] / detmech::W[k];
            }
            reaction_rates(rho, Y, T0, C0, om0);
            Real om_max = 0.0_rt;
            for (int k = 0; k < NS; ++k) {
                om_max = std::max(om_max, std::abs(om0[k]));
            }
            out[NS + 6] = om_max;
            Real omj[NS], dTv[NS], um[NS], jj[NS][NS];
            Real cvv = 0.0_rt;
            reaction_rates_jac(T0, C0, omj, jj, dTv, um, cvv);
            Real omj_max = 0.0_rt;
            for (int k = 0; k < NS; ++k) {
                omj_max = std::max(omj_max, std::abs(omj[k]));
            }
            out[NS + 7] = omj_max;
            Real T1 = T0;
            Real dy = 0.0_rt;
            Real dT = 0.0_rt;
            out[NS + 8] = 0.0_rt;
            out[NS + 9] = 0.0_rt;
            out[NS + 10] = 0.0_rt;
            const bool implicit_ok = implicit_chemistry_step(
              C0, rho * e0, T0, 1.0e-8_rt, 0.05_rt, 0.03_rt,
              250.0_rt, 6500.0_rt, C1, T1, dy, dT, out + NS + 8);
            out[NS] = implicit_ok ? 1.0_rt : 0.0_rt;
            out[NS + 1] = dy;
            out[NS + 2] = dT;

            if (!implicit_ok) {
                for (int k = 0; k < NS; ++k) {
                    C1[k] = C0[k];
                }
                T1 = T0;
            }
            Real sum_y = 0.0_rt;
            for (int k = 0; k < NS; ++k) {
                Y[k] = C1[k] * detmech::W[k] / rho;
                sum_y += Y[k];
            }
            for (int k = 0; k < NS; ++k) {
                Y[k] /= sum_y;
            }
            Real T2 = T1;
            Real hnext = 1.0e-8_rt;
            Real activity = 0.0_rt;
            chemistry_step(
              rho, e0, 4.64e-6_rt, 650.0_rt, 1.0e-8_rt,
              0.05_rt, 0.03_rt, Y, T2, hnext, activity);
            for (int k = 0; k < NS; ++k) {
                out[k] = Y[k];
            }
            out[NS + 3] = T2;
            out[NS + 4] = rho;
            out[NS + 5] = e0;
        });
        amrex::Gpu::streamSynchronize();

        std::vector<amrex::Real> host(nout);
        amrex::Gpu::copyAsync(
          amrex::Gpu::deviceToHost, device.begin(), device.end(),
          host.begin());
        amrex::Gpu::streamSynchronize();
        std::cout << std::setprecision(10);
        std::cout << "implicit_ok=" << host[det::NS]
                  << " dy=" << host[det::NS + 1]
                  << " dT=" << host[det::NS + 2]
                  << " omega_max=" << host[det::NS + 6]
                  << " omj_max=" << host[det::NS + 7]
                  << " dbg_it=" << host[det::NS + 8]
                  << " dbg_res=" << host[det::NS + 9]
                  << " dbg_fail=" << host[det::NS + 10] << '\n';
        for (int k = 0; k < det::NS; ++k) {
            std::cout << detmech::species_names[k] << "=" << host[k] << '\n';
        }
        std::cout << "T=" << host[det::NS + 3]
                  << " rho=" << host[det::NS + 4]
                  << " e0=" << host[det::NS + 5] << '\n';
    }
    amrex::Finalize();
}
