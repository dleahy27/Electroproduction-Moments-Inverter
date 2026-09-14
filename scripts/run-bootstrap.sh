#!/usr/bin/env bash

# Might as well call the build/emi in case someone forgets to alias
# Built in redundancy essentially
project_dir="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

# Command line option for toys etc, mainly just use default
options=()
[[ -n "${TOYS:-}" ]] && options+=(--toys "${TOYS}")
[[ -n "${STARTS_PER_TOY:-}" ]] && options+=(--starts-per-toy "${STARTS_PER_TOY}")
[[ -n "${WORKERS:-}" ]] && options+=(--workers "${WORKERS}")

# Here you can change what bin(s) and file you want to look at, also the value for epsilon
"${project_dir}/build/emi" bootstrap "${options[@]}" \
  --electro \
  --input "${project_dir}/InputFiles/Experiment/e_rho_moments.root" \
  --tree expMoments \
  --bin "${1:-0}" \
  --epsilon 0.8 \
  --output "${project_dir}/OutputFiles/e_rho_bootstrap_${1:-0}.root"
