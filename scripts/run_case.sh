#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_file="${1:-${project_dir}/examples/h2o2_2d_demo/inputs}"
run_dir="${DETONATION_RUN_DIR:-${project_dir}/runs/demo}"

if [[ "${input_file}" != /* ]]; then
  input_file="${project_dir}/${input_file}"
fi
if [[ "${run_dir}" != /* ]]; then
  run_dir="${project_dir}/${run_dir}"
fi
if [[ ! -f "${input_file}" ]]; then
  echo "Input file not found: ${input_file}" >&2
  exit 1
fi

mkdir -p "${run_dir}"
shift || true

mpi_cmd=(mpirun)
if [[ "$(id -u)" -eq 0 ]]; then
  mpi_cmd+=(--allow-run-as-root)
fi
mpi_cmd+=(-n "${NP:-1}")

cd "${run_dir}"
"${mpi_cmd[@]}" \
  "${project_dir}/build/detonation2d" "${input_file}" "$@" \
  2>&1 | tee -a run.log
