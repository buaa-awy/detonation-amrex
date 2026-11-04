#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
amrex_prefix="${AMREX_PREFIX:-${HOME}/amrex-install-2d-cuda}"
cuda_home="${CUDA_HOME:-}"

if [[ ! -f "${amrex_prefix}/lib/cmake/AMReX/AMReXConfig.cmake" ]]; then
  echo "AMReX was not found at ${amrex_prefix}" >&2
  echo "Set AMREX_PREFIX to a 2-D AMReX CUDA installation." >&2
  exit 1
fi

if [[ -z "${cuda_home}" ]]; then
  for candidate in /usr/local/cuda /usr/local/cuda-12.9; do
    if [[ -x "${candidate}/bin/nvcc" ]]; then
      cuda_home="${candidate}"
      break
    fi
  done
fi
if [[ -z "${cuda_home}" || ! -x "${cuda_home}/bin/nvcc" ]]; then
  echo "nvcc was not found. Set CUDA_HOME to the CUDA toolkit directory." >&2
  exit 1
fi

export PATH="${cuda_home}/bin:${PATH}"
cuda_arch="${CUDA_ARCH:-89}"
jobs="${BUILD_JOBS:-2}"
cxx_compiler="${CXX_COMPILER:-g++}"

cmake -S "${project_dir}" -B "${project_dir}/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="${cxx_compiler}" \
  -DCMAKE_CUDA_COMPILER="${cuda_home}/bin/nvcc" \
  -DDETONATION_CUDA_ARCHITECTURES="${cuda_arch}" \
  -DCMAKE_PREFIX_PATH="${amrex_prefix}"

cmake --build "${project_dir}/build" --parallel "${jobs}"
