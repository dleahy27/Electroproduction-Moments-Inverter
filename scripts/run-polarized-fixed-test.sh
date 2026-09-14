#!/usr/bin/env bash

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
emi="${project_dir}/build/emi"

if [[ ! -x "${emi}" ]]; then
  echo "${emi} does not exist; run scripts/build.sh first" >&2
  exit 1
fi

epsilon="${EPSILON:-0.8}"
starts="${STARTS:-10000}"
workers="${WORKERS:-1}"
seed="${SEED:-12345}"
gauge="${K_GAUGE:-}"

mkdir -p "${project_dir}/InputFiles/Generated" "${project_dir}/OutputFiles"

for mode in initial recoil both; do
  input="${project_dir}/InputFiles/Generated/fixed_${mode}_moments.root"
  output="${project_dir}/OutputFiles/fixed_${mode}_fit.root"
  gauge_options=()
  if [[ -n "${gauge}" && "${mode}" != "both" ]]; then
    gauge_options=(--k-gauge "${gauge}")
  fi

  echo "Generating ${mode}-polarization closure input"
  "${emi}" generate-example \
    --polarization "${mode}" \
    --epsilon "${epsilon}" \
    --output "${input}" \
    --quiet \
    "${gauge_options[@]}"

  fit_options=()
  if [[ "${HESSE:-0}" != "1" ]]; then
    fit_options=(--no-hesse)
  fi

  echo "Fitting ${mode}-polarization closure input"
  "${emi}" fit \
    --electro \
    --input "${input}" \
    --tree genMoments \
    --output "${output}" \
    --polarization "${mode}" \
    --epsilon "${epsilon}" \
    --starts "${starts}" \
    --workers "${workers}" \
    --seed "${seed}" \
    --quiet \
    "${gauge_options[@]}" \
    "${fit_options[@]}"
done
