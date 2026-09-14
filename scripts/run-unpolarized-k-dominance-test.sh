#!/usr/bin/env bash

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
emi="${project_dir}/build/emi"

if [[ ! -x "${emi}" ]]; then
  echo "${emi} does not exist; run scripts/build.sh first" >&2
  exit 1
fi

settings="${SETTINGS:-${project_dir}/app/UserSettings.h}"
epsilon="${EPSILON:-}"
starts="${STARTS:-1000}"
workers="${WORKERS:-8}"
seed="${SEED:-0}"
repeats="${REPEATS:-200}"
suppression_text="${SUPPRESSIONS:-}"

if [[ ! "${repeats}" =~ ^[1-9][0-9]*$ ]]; then
  echo "REPEATS must be a positive integer, got '${repeats}'" >&2
  exit 1
fi

if [[ ! -f "${settings}" ]]; then
  echo "Settings file does not exist: ${settings}" >&2
  exit 1
fi

if [[ -n "${suppression_text}" ]]; then
  IFS=', ' read -r -a suppressions <<< "${suppression_text}"
else
  suppressions=("")
fi

batch="N${repeats}"
if [[ -n "${suppression_text}" ]]; then
  batch="${batch}_M${#suppressions[@]}"
fi
input_dir="${project_dir}/InputFiles/Generated/photo_k/${batch}"
output_dir="${project_dir}/OutputFiles/photo_k/${batch}"
mkdir -p "${input_dir}" "${output_dir}"

fit_options=()
if [[ "${HESSE:-0}" != "1" ]]; then
  fit_options=(--no-hesse)
fi

for m in "${!suppressions[@]}"; do
  suppression="${suppressions[m]}"
  for ((n = 0; n < repeats; ++n)); do
    run_seed=$((seed + m * repeats + n))
    stem="fixed_N${n}"
    if [[ -n "${suppression}" ]]; then
      stem="${stem}_M${m}"
    fi
    input="${input_dir}/${stem}_moments.root"
    output="${output_dir}/${stem}_unpolarized_fit.root"

    generation_options=(
      --settings "${settings}"
      --output "${input}"
      --seed "${run_seed}"
      --quiet
    )
    if [[ -n "${epsilon}" ]]; then
      generation_options+=(--epsilon "${epsilon}")
    fi
    if [[ -n "${suppression}" ]]; then
      generation_options+=(--suppression "${suppression}")
    fi

    echo "Generating configured fixed point: N=${n}, seed=${run_seed}${suppression:+, suppression=${suppression}}"
    "${emi}" generate-fixed "${generation_options[@]}"

    echo "Fitting its unpolarized moments: N=${n}"
    fit_runtime_options=()
    if [[ -n "${epsilon}" ]]; then
      fit_runtime_options+=(--epsilon "${epsilon}")
    fi
    "${emi}" fit \
      --settings "${settings}" \
      --input "${input}" \
      --tree genMoments \
      --output "${output}" \
      --polarization none \
      --starts "${starts}" \
      --workers "${workers}" \
      --seed "${run_seed}" \
      --quiet \
      "${fit_runtime_options[@]}" \
      "${fit_options[@]}"
  done
done

echo "Wrote generated files under ${input_dir}"
echo "Wrote fitted files under ${output_dir}"
echo "For a single result, compare it with AnalysisScripts/polarized_fixed_amp_analysis.ipynb."
