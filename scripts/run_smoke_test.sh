#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
stamp="$(date +%Y%m%d-%H%M%S)"
run_dir="${DETONATION_SMOKE_DIR:-${project_dir}/runs/smoke-${stamp}}"

DETONATION_RUN_DIR="${run_dir}" \
  "${project_dir}/scripts/run_case.sh" \
  "${project_dir}/examples/h2o2_2d_demo/inputs-smoke"

python3 "${project_dir}/scripts/check_smoke.py" \
  "${run_dir}/status.json"
echo "Smoke run directory: ${run_dir}"
