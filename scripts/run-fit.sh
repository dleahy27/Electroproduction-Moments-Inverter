#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

options=()
[[ -n "${STARTS:-}" ]] && options+=(--starts "${STARTS}")
[[ -n "${WORKERS:-}" ]] && options+=(--workers "${WORKERS}")

"${project_dir}/build/emi" fit "${options[@]}" \
  --input "${project_dir}/InputFiles/Experiment/e_rho_moments.root" \
  --tree expMoments \
  --bin "${1:-0}" \
  --epsilon 0.8 \
  --output "${project_dir}/OutputFiles/e_rho_fit_${1:-0}.root"
