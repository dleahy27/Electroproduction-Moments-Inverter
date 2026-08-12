#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

options=()
[[ -n "${TOYS:-}" ]] && options+=(--toys "${TOYS}")
[[ -n "${STARTS_PER_TOY:-}" ]] && options+=(--starts-per-toy "${STARTS_PER_TOY}")
[[ -n "${WORKERS:-}" ]] && options+=(--workers "${WORKERS}")

"${project_dir}/build/emi" bootstrap "${options[@]}" \
  --input "${project_dir}/InputFiles/Experiment/e_rho_moments.root" \
  --tree expMoments \
  --bin "${1:-0}" \
  --epsilon 0.8 \
  --output "${project_dir}/OutputFiles/e_rho_bootstrap_${1:-0}.root"
