#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bin="${1:-0}"
starts="${2:-100}"
common=(
  --input "${project_dir}/InputFiles/Experiment/e_rho_moments.root"
  --tree expMoments
  --bin "${bin}"
  --epsilon 0.8
  --starts "${starts}"
  --workers 1
  --seed 12345
  --no-hesse
)

echo "Analytical gradients"
time "${project_dir}/build/emi" fit "${common[@]}" \
  --output "${project_dir}/OutputFiles/gradient_analytic_${bin}.root"

echo "Numerical gradients"
time "${project_dir}/build/emi" fit "${common[@]}" \
  --numerical-gradients \
  --output "${project_dir}/OutputFiles/gradient_numerical_${bin}.root"
