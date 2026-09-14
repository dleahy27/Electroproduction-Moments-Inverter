#!/usr/bin/env bash

# Define all the parameters that are needed to run scripts automatically
RHO_BINS=3
OMEGA_BINS=2
GLUEX_BINS=17

E_RHO_EPS=0.8
MU_RHO_EPS=0.9
OMEGA_EPS=0.96

# cml args (default should be fine for most part is what I used for results)
options=()
[[ -n "${TOYS:-}" ]] && options+=(--toys "${TOYS}")
[[ -n "${STARTS_PER_TOY:-}" ]] && options+=(--starts-per-toy "${STARTS_PER_TOY}")
[[ -n "${WORKERS:-}" ]] && options+=(--workers "${WORKERS}")

for i in $(seq 0 $RHO_BINS); do
  "${project_dir}/build/emi" bootstrap "${options[@]}" \
    --input "${project_dir}/InputFiles/Experiment/e_rho_moments.root" \
    --tree expMoments \
    --bin "${i}" \
    --epsilon E_RHO_EPS \
    --output "${project_dir}/OutputFiles/e_rho_bootstrap_${i}.root"

    "${project_dir}/build/emi" bootstrap "${options[@]}" \
        --input "${project_dir}/InputFiles/Experiment/mu_rho_moments.root" \
        --tree expMoments \
        --bin "${i}" \
        --epsilon MU_RHO_EPS \
        --output "${project_dir}/OutputFiles/mu_rho_bootstrap_${i}.root"
done

for i in $(seq 0 $OMEGA_BINS); do
  "${project_dir}/build/emi" bootstrap "${options[@]}" \
    --input "${project_dir}/InputFiles/Experiment/e_omega_moments.root" \
    --tree expMoments \
    --bin "${i}" \
    # Eps is the same for e_omega as for e_rho
    --epsilon E_RHO_EPS \
    --output "${project_dir}/OutputFiles/e_omega_bootstrap_${i}.root"

    "${project_dir}/build/emi" bootstrap "${options[@]}" \
        --input "${project_dir}/InputFiles/Experiment/mu_omega_moments.root" \
        --tree expMoments \
        --bin "${i}" \
        --epsilon OMEGA_EPS \
        --output "${project_dir}/OutputFiles/mu_omega_bootstrap_${i}.root"

    # Phi data, is in supplementary notes but has not been fully published, just in a thesis
#    "${project_dir}/build/emi" bootstrap "${options[@]}" \
#            --input "${project_dir}/InputFiles/Experiment/e_phi_moments.root" \
#            --tree expMoments \
#            --bin "${i}" \
#            --epsilon E_RHO_EPS \
#            --output "${project_dir}/OutputFiles/e_phi_bootstrap_${i}.root"
done

# GlueX photoprod bootstrap plots, given by Derek but not focus
# Need to double check if refactor can handle these... zeroed out moments should be ignored
#for i in $(seq 0 ${GLUEX_BINS}); do
#  "${project_dir}/build/emi" bootstrap "${options[@]}" \
#    --input "${project_dir}/InputFiles/Experiment/gluex_moments.root" \
#    --tree expMoments \
#    --bin "${i}" \
#    --epsilon 0 \
#    --output "${project_dir}/OutputFiles/gluex_bootstrap_${i}.root"
#done