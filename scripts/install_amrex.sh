#!/usr/bin/env bash
set -euo pipefail

amrex_version="${AMREX_VERSION:-26.10}"
amrex_commit="${AMREX_COMMIT:-015c14b65715b277fbc54706ebe7455b196d2254}"
amrex_repo="${AMREX_REPO:-https://github.com/AMReX-Codes/amrex.git}"
amrex_source="${AMREX_SOURCE:-${HOME}/amrex}"
amrex_build="${AMREX_BUILD:-${HOME}/amrex-build-2d-cuda}"
amrex_prefix="${AMREX_PREFIX:-${HOME}/amrex-install-2d-cuda}"
cuda_home="${CUDA_HOME:-}"
cuda_arch="${CUDA_ARCH:-89}"
jobs="${BUILD_JOBS:-4}"

if command -v git >/dev/null 2>&1; then :; else
  echo "git was not found." >&2
  exit 1
fi
if command -v cmake >/dev/null 2>&1; then :; else
  echo "CMake was not found." >&2
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
  echo "nvcc was not found. Set CUDA_HOME to the CUDA Toolkit directory." >&2
  exit 1
fi

if [[ ! -d "${amrex_source}/.git" ]]; then
  if [[ -e "${amrex_source}" ]]; then
    echo "${amrex_source} exists but is not a Git checkout." >&2
    echo "Choose an empty AMREX_SOURCE path." >&2
    exit 1
  fi
  git clone --branch "${amrex_version}" --depth 1 \
    "${amrex_repo}" "${amrex_source}"
else
  current="$(git -C "${amrex_source}" describe --tags --always 2>/dev/null || true)"
  if [[ "${current}" != "${amrex_version}" ]]; then
    echo "AMReX at ${amrex_source} is '${current}', expected '${amrex_version}'." >&2
    echo "Set AMREX_SOURCE to a fresh path for a reproducible installation." >&2
    exit 1
  fi
  if ! git -C "${amrex_source}" diff --quiet || \
     ! git -C "${amrex_source}" diff --cached --quiet; then
    echo "AMReX source has tracked modifications: ${amrex_source}" >&2
    exit 1
  fi
fi

actual_commit="$(git -C "${amrex_source}" rev-parse HEAD)"
if [[ "${actual_commit}" != "${amrex_commit}" ]]; then
  echo "AMReX commit is ${actual_commit}, expected ${amrex_commit}." >&2
  exit 1
fi

cmake -S "${amrex_source}" -B "${amrex_build}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAMReX_SPACEDIM=2 \
  -DAMReX_GPU_BACKEND=CUDA \
  -DCMAKE_CUDA_ARCHITECTURES="${cuda_arch}" \
  -DCMAKE_CUDA_COMPILER="${cuda_home}/bin/nvcc" \
  -DAMReX_MPI=ON \
  -DAMReX_OMP=OFF \
  -DAMReX_EB=OFF \
  -DAMReX_PARTICLES=OFF \
  -DAMReX_PLOTFILE_TOOLS=ON \
  -DCMAKE_INSTALL_PREFIX="${amrex_prefix}"

cmake --build "${amrex_build}" --parallel "${jobs}"
cmake --install "${amrex_build}"

fextract="${amrex_prefix}/bin/amrex_fextract"
if [[ ! -x "${fextract}" ]]; then
  echo "AMReX installed, but ${fextract} was not generated." >&2
  exit 1
fi

cat <<EOF
AMReX ${amrex_version} installed at:
  ${amrex_prefix}

For this shell:
  export AMREX_PREFIX=${amrex_prefix}
  export CUDA_HOME=${cuda_home}
  export CUDA_ARCH=${cuda_arch}
  export FEXTRACT=${fextract}
EOF
