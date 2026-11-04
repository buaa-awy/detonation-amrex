#include <AMReX.H>
#include <AMReX_Gpu.H>

#include "Detonation.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::Print() << "AMReX 2D CUDA device: "
                       << amrex::Gpu::Device::deviceName() << '\n';
        Detonation solver;
        solver.InitData();
        solver.Evolve();
    }
    amrex::Finalize();
    return 0;
}
