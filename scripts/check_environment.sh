#!/usr/bin/env bash
set -u

status=0
check_command() {
  local name="$1"
  if command -v "${name}" >/dev/null 2>&1; then
    printf 'OK      %-12s %s\n' "${name}" "$(command -v "${name}")"
  else
    printf 'MISSING %-12s\n' "${name}"
    status=1
  fi
}

check_command git
check_command g++
check_command cmake
check_command mpirun
check_command python3

cuda_home="${CUDA_HOME:-/usr/local/cuda}"
if [[ -x "${cuda_home}/bin/nvcc" ]]; then
  printf 'OK      %-12s %s\n' nvcc "${cuda_home}/bin/nvcc"
else
  printf 'MISSING %-12s %s\n' nvcc "${cuda_home}/bin/nvcc"
  status=1
fi

amrex_prefix="${AMREX_PREFIX:-${HOME}/amrex-install-2d-cuda}"
if [[ -f "${amrex_prefix}/lib/cmake/AMReX/AMReXConfig.cmake" ]]; then
  printf 'OK      %-12s %s\n' AMReX "${amrex_prefix}"
else
  printf 'MISSING %-12s %s\n' AMReX "${amrex_prefix}"
  status=1
fi

fextract="${FEXTRACT:-${amrex_prefix}/bin/amrex_fextract}"
if [[ -x "${fextract}" ]]; then
  printf 'OK      %-12s %s\n' fextract "${fextract}"
else
  printf 'WARNING %-12s %s\n' fextract "not found"
fi

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,driver_version,memory.total \
    --format=csv,noheader
fi

exit "${status}"
