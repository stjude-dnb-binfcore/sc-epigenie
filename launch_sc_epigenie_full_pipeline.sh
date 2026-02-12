#!/bin/bash
###############################################################################
# Script:        launch_sc_epigenie_full_pipeline.sh
# Purpose:       Launch the full sc-epigenie pipeline on LSF as a dependency-
#                aware job chain, with fully OPTIONAL steps and adaptive waits.
#
# Summary:
#   - Submits the end-to-end sc-epigenie workflow (A→I) to LSF.
#   - Every step can be toggled on/off independently.
#   - Dependencies automatically wire up only between ENABLED steps.
#   - CellRanger can run in parallel with FastQC or wait for it, based on a flag.
#   - Email notifications sent on completion (and optionally on start) of each job.
#
# Key Features:
#   • Full workflow coverage with modular, optional steps:
#       A: FastQC
#       B: CellRanger
#       C: Upstream analysis
#       D: Integrative analysis
#       E: Cluster cell calling
#       F: Integration with scRNA-seq
#       G: Trajectories (Monocle)
#       H: Motif footprint analysis
#       I: R/Shiny app packaging
#   • Optional per-step execution via RUN_* flags.
#   • Adaptive dependencies: each enabled step depends on the last enabled predecessor.
#   • LSF email notifications with a single NOTIFY_EMAIL value for all jobs.
#
# Usage:
#   1) Edit the configuration section near the top:
#        - NOTIFY_EMAIL="user.name@stjude.org"
#        - NOTIFY_ON_START=0|1
#        - RUN_FASTQC / RUN_CELLRANGER / ... / RUN_RSHINY = 0|1
#   2) Ensure each module directory has its LSF script:
#        analyses/<module>/lsf-script.txt
#        (CellRanger may use submit-multiple-jobs.sh to fan out internal jobs.)
#   3) Run from anywhere:
#        bash launch_sc_epigenie_full_pipeline.sh
#
# Notifications:
#   - Completion emails: enabled by default via `-N -u "$NOTIFY_EMAIL"`.
#   - Start emails: set NOTIFY_ON_START=1 to add `-B` to all bsubs.
#   - Cluster defaults like LSB_MAILTO may be honored but are overridden by -u.
#
# Dependencies Semantics:
#   - The chain uses `done(<JOBID>)` so downstream steps start only if upstream
#     succeeded. If you want downstream steps to run regardless of status,
#     switch to `ended(<JOBID>)` where appropriate.
#   - If you need a step to wait on multiple predecessors, combine expressions:
#       -w "done(<JOB_A>) && done(<JOB_B>)"
#
# Logs:
#   - Each step writes to <module>/job.out and <module>/job.err.
#   - LSF submission output is parsed to capture Job IDs for a final summary.
#
# Requirements:
#   - LSF (bsub) available on PATH.
#   - Module directories exist with their respective LSF scripts.
#
# Maintainer:    Antonia Chroni (DNB Bioinformatics, St. Jude Children's Research Hospital)
# Last updated:  2026-02-09
###############################################################################

set -e
set -o pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
PROJECT_DIR="$(pwd)"

# ------------------------------------------------------------------------------
# Notifications
# ------------------------------------------------------------------------------
# All jobs will notify on completion (-N) to this email.
# Set NOTIFY_ON_START=1 to also get a mail when a job starts (-B).
NOTIFY_EMAIL="user.name@stjude.org"
NOTIFY_ON_START=1   # 1 = include -B (mail on start); 0 = Skip start notifications

# Compose common bsub notification flags
BSUB_NOTIFY_FLAGS=(-N -u "${NOTIFY_EMAIL}")
if [[ "${NOTIFY_ON_START}" -eq 1 ]]; then
  BSUB_NOTIFY_FLAGS=(-B "${BSUB_NOTIFY_FLAGS[@]}")
fi

# ------------------------------------------------------------------------------
# Feature toggles — make EVERY step optional
# 1 = run step; 0 = skip step
# ------------------------------------------------------------------------------
RUN_FASTQC=1                # A: FastQC: `fastqc-analysis`
RUN_CELLRANGER=1            # B: CellRanger: `cellranger-analysis`
RUN_UPSTREAM=1              # C: Upstream: `upstream-analysis`
RUN_INTEGRATIVE=1           # D: Integrative: `integrative-analysis`
RUN_CLUSTER=1               # E: Cluster cell calling: `cluster-cell-calling`
RUN_INTEGRATE_SCRNA=1       # F: Integration with scRNA-seq: `integration-with-scrna-seq-data`
RUN_TRAJECTORIES=1          # G: Trajectories (Monocle): `building-trajectories-monocle`
RUN_MOTIF=1                 # H: Motif footprint analysis: `motif-footprint-analysis`
RUN_RSHINY=1                # I: R/Shiny app: `rshiny-app`

# ------------------------------------------------------------------------------
# Heading
# ------------------------------------------------------------------------------
echo "============================================================"
echo " Launching sc-epigenie pipeline (LSF job chain) — optional steps"
echo " Working directory: ${PROJECT_DIR}"
echo "============================================================"

# ------------------------------------------------------------------------------
# Helper: extract LSF job ID from 'Job <12345> is submitted ...'
# ------------------------------------------------------------------------------
extract_job_id() { awk '{print $2}' | sed 's/[<>]//g'; }

# ------------------------------------------------------------------------------
# Helper: bsub wrapper to reduce repetition
#   args:
#     1 = workdir
#     2 = stdin script path
#     3 = dependency expression (may be empty)
#     4 = label for logs
#   echoes job id
# ------------------------------------------------------------------------------
submit_job() {
  local _cwd="$1"
  local _in="$2"
  local _dep="${3:-}"
  local _label="${4:-job}"

  mkdir -p "${_cwd}"
  local out
  if ! out=$(bsub -cwd "${_cwd}" \
                  ${_dep:+-w "${_dep}"} \
                  "${BSUB_NOTIFY_FLAGS[@]}" \
                  < "${_in}" 2>&1); then
    echo "ERROR: bsub failed for ${_label} in ${_cwd}" >&2
    echo "${out}" >&2
    exit 1
  fi
  local jid
  jid="$(echo "${out}" | extract_job_id)"
  if [[ -z "${jid}" ]]; then
    echo "ERROR: Could not parse job ID for ${_label}. bsub output:" >&2
    echo "${out}" >&2
    exit 1
  fi
  echo "${jid}"
}

# ------------------------------------------------------------------------------
# Helper: bsub wrapper for CellRanger's submit-multiple-jobs.sh, which may print multiple lines of output.

# submit_job WORKDIR LAUNCHER_SCRIPT UPSTREAM_DEP JOBNAME
# Example:
#   JOB_B=$(submit_job "${B_DIR}" "${B_DIR}/submit-multiple-jobs.sh" "${B_DEP}" "CellRanger")
#
# Contract:
#   - LAUNCHER_SCRIPT runs LOCALLY (not as an LSF job).
#   - It is responsible for submitting all internal LSF jobs and their own waiter.
#   - It must print ONLY the final waiter job ID (numeric) to stdout.
submit_job_cellranger() {
  local workdir="$1"           # e.g., ${B_DIR}
  local launcher="$2"          # e.g., ${B_DIR}/submit-multiple-jobs.sh
  local upstream_dep="$3"      # e.g., "done(${JOB_A})" or ""
  local jobname="$4"           # logical prefix, e.g., "CellRanger"

  (
    set -euo pipefail
    cd "$workdir"

    # Run launcher LOCALLY. Capture ALL output.
    local raw
    raw="$("${launcher}" "${upstream_dep}" "${jobname}" 2>&1)"

    # Extract the LAST numeric token from output (the waiter job ID)
    local waiter_id
    waiter_id="$(echo "${raw}" | grep -Eo 'Job <[0-9]+>' | grep -Eo '[0-9]+' | tail -n 1)"

    # Validate
    if [[ -z "${waiter_id}" || ! "${waiter_id}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: ${launcher} did not emit a numeric job ID." >&2
      echo "       Raw launcher output was:" >&2
      echo "-----------------------------------------------" >&2
      echo "${raw}" >&2
      echo "-----------------------------------------------" >&2
      exit 1
    fi

    # Return ONLY the numeric waiter job ID
    echo "${waiter_id}"
  )
}

# ------------------------------------------------------------------------------
# Module directories
# ------------------------------------------------------------------------------
A_DIR="${PROJECT_DIR}/analyses/fastqc-analysis"
B_DIR="${PROJECT_DIR}/analyses/cellranger-analysis"
C_DIR="${PROJECT_DIR}/analyses/upstream-analysis"
D_DIR="${PROJECT_DIR}/analyses/integrative-analysis"
E_DIR="${PROJECT_DIR}/analyses/cluster-cell-calling"
F_DIR="${PROJECT_DIR}/analyses/integration-with-scrna-seq-data"
G_DIR="${PROJECT_DIR}/analyses/building-trajectories-monocle"
H_DIR="${PROJECT_DIR}/analyses/motif-footprint-analysis"
I_DIR="${PROJECT_DIR}/analyses/rshiny-app"

# ------------------------------------------------------------------------------
# Dynamic progress counter (counts only enabled steps)
# ------------------------------------------------------------------------------
count_enabled() {
  local n=0
  (( RUN_FASTQC ))           && ((n++))
  (( RUN_CELLRANGER ))       && ((n++))
  (( RUN_UPSTREAM ))         && ((n++))
  (( RUN_INTEGRATIVE ))      && ((n++))
  (( RUN_CLUSTER ))          && ((n++))
  (( RUN_INTEGRATE_SCRNA ))  && ((n++))
  (( RUN_TRAJECTORIES ))     && ((n++))
  (( RUN_MOTIF ))            && ((n++))
  (( RUN_RSHINY ))           && ((n++))
  echo "${n}"
}
TOTAL_STEPS="$(count_enabled)"
STEP=0

bump_step() { STEP=$((STEP+1)); }

# ------------------------------------------------------------------------------
# Submission chain
# ------------------------------------------------------------------------------
LAST_JOB=""
JOB_A=""; JOB_B=""; JOB_C=""; JOB_D=""; JOB_E=""; JOB_F=""; JOB_G=""; JOB_H=""; JOB_I=""

# A) FastQC
if (( RUN_FASTQC )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting FastQC (A)..."
  JOB_A=$(submit_job "${A_DIR}" "${A_DIR}/lsf-script.txt" "" "FastQC")
  echo "  A(FastQC) = ${JOB_A}"
else
  echo "[–/–] FastQC (A): SKIPPED"
fi


# B) CellRanger 
if (( RUN_CELLRANGER )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting CellRanger (B)..."
  B_DEP=""

  # submit-multiple-jobs.sh 
  JOB_B=$(submit_job_cellranger "${B_DIR}" "${B_DIR}/submit-multiple-jobs.sh" "${B_DEP}" "CellRanger")
  # Email me when the CellRanger job (JOB_B) is submitted.
  echo "  B(CellRanger) = ${JOB_B} ${B_DEP:+(dep: ${B_DEP})}" | mail -s "CellRanger submitted" "${NOTIFY_EMAIL}"
else
  echo "[–/–] CellRanger (B): SKIPPED"
fi


# Maintain a simple linear chain among *enabled* steps after B:
# The next step depends on the most recent ENABLED job. If B ran, chain from B; else if A ran, chain from A; else no dep.
if   (( RUN_CELLRANGER )); then LAST_JOB="${JOB_B}"
elif (( RUN_FASTQC    )); then LAST_JOB="${JOB_A}"
else                            LAST_JOB=""
fi

# C) Upstream
if (( RUN_UPSTREAM )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Upstream (C)..."
  C_DEP=""
  [[ -n "${LAST_JOB}" ]] && C_DEP="done(${LAST_JOB})"
  JOB_C=$(submit_job "${C_DIR}" "${C_DIR}/lsf-script.txt" "${C_DEP}" "Upstream")
  echo "  C(Upstream) = ${JOB_C} ${C_DEP:+(dep: ${C_DEP})}"
  LAST_JOB="${JOB_C}"
else
  echo "[–/–] Upstream (C): SKIPPED"
fi

# D) Integrative
if (( RUN_INTEGRATIVE )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Integrative (D)..."
  D_DEP=""
  [[ -n "${LAST_JOB}" ]] && D_DEP="done(${LAST_JOB})"
  JOB_D=$(submit_job "${D_DIR}" "${D_DIR}/lsf-script.txt" "${D_DEP}" "Integrative")
  echo "  D(Integrative) = ${JOB_D} ${D_DEP:+(dep: ${D_DEP})}"
  LAST_JOB="${JOB_D}"
else
  echo "[–/–] Integrative (D): SKIPPED"
fi

# E) Cluster cell calling
if (( RUN_CLUSTER )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Cluster cell calling (E)..."
  E_DEP=""
  [[ -n "${LAST_JOB}" ]] && E_DEP="done(${LAST_JOB})"
  JOB_E=$(submit_job "${E_DIR}" "${E_DIR}/lsf-script.txt" "${E_DEP}" "Cluster cell calling")
  echo "  E(Cluster) = ${JOB_E} ${E_DEP:+(dep: ${E_DEP})}"
  LAST_JOB="${JOB_E}"
else
  echo "[–/–] Cluster cell calling (E): SKIPPED"
fi

# F) Integration with scRNA-seq
if (( RUN_INTEGRATE_SCRNA )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Integration with scRNA-seq (F)..."
  F_DEP=""
  [[ -n "${LAST_JOB}" ]] && F_DEP="done(${LAST_JOB})"
  JOB_F=$(submit_job "${F_DIR}" "${F_DIR}/lsf-script.txt" "${F_DEP}" "Integration with scRNA-seq")
  echo "  F(Integration scRNA) = ${JOB_F} ${F_DEP:+(dep: ${F_DEP})}"
  LAST_JOB="${JOB_F}"
else
  echo "[–/–] Integration with scRNA-seq (F): SKIPPED"
fi

# G) Trajectories (Monocle)
if (( RUN_TRAJECTORIES )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Trajectories (Monocle) (G)..."
  G_DEP=""
  [[ -n "${LAST_JOB}" ]] && G_DEP="done(${LAST_JOB})"
  JOB_G=$(submit_job "${G_DIR}" "${G_DIR}/lsf-script.txt" "${G_DEP}" "Trajectories (Monocle)")
  echo "  G(Monocle) = ${JOB_G} ${G_DEP:+(dep: ${G_DEP})}"
  LAST_JOB="${JOB_G}"
else
  echo "[–/–] Trajectories (Monocle) (G): SKIPPED"
fi

# H) Motif footprint analysis
if (( RUN_MOTIF )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting Motif footprint analysis (H)..."
  H_DEP=""
  [[ -n "${LAST_JOB}" ]] && H_DEP="done(${LAST_JOB})"
  JOB_H=$(submit_job "${H_DIR}" "${H_DIR}/lsf-script.txt" "${H_DEP}" "Motif footprint")
  echo "  H(Motif) = ${JOB_H} ${H_DEP:+(dep: ${H_DEP})}"
  LAST_JOB="${JOB_H}"
else
  echo "[–/–] Motif footprint (H): SKIPPED"
fi

# I) R/Shiny app
if (( RUN_RSHINY )); then
  bump_step
  echo "[${STEP}/${TOTAL_STEPS}] Submitting R/Shiny app (I)..."
  I_DEP=""
  [[ -n "${LAST_JOB}" ]] && I_DEP="done(${LAST_JOB})"
  JOB_I=$(submit_job "${I_DIR}" "${I_DIR}/lsf-script.txt" "${I_DEP}" "R/Shiny app")
  echo "  I(R/Shiny) = ${JOB_I} ${I_DEP:+(dep: ${I_DEP})}"
  LAST_JOB="${JOB_I}"
else
  echo "[–/–] R/Shiny app (I): SKIPPED"
fi

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo "============================================================"
if (( TOTAL_STEPS == 0 )); then
  echo " No steps were enabled — nothing submitted."
else
  echo " All enabled steps submitted."
fi
printf "  A(FastQC)   : %s\n"   "${RUN_FASTQC:+${JOB_A:-SKIPPED}}"
printf "  B(CellRanger): %s\n"   "${RUN_CELLRANGER:+${JOB_B:-SKIPPED}}"
printf "  C(Upstream) : %s\n"    "${RUN_UPSTREAM:+${JOB_C:-SKIPPED}}"
printf "  D(Integrative): %s\n"  "${RUN_INTEGRATIVE:+${JOB_D:-SKIPPED}}"
printf "  E(Cluster)  : %s\n"    "${RUN_CLUSTER:+${JOB_E:-SKIPPED}}"
printf "  F(Integr-scRNA): %s\n" "${RUN_INTEGRATE_SCRNA:+${JOB_F:-SKIPPED}}"
printf "  G(Monocle)  : %s\n"    "${RUN_TRAJECTORIES:+${JOB_G:-SKIPPED}}"
printf "  H(Motif)    : %s\n"    "${RUN_MOTIF:+${JOB_H:-SKIPPED}}"
printf "  I(R/Shiny)  : %s\n"    "${RUN_RSHINY:+${JOB_I:-SKIPPED}}"
echo " Final job in chain: ${LAST_JOB:-NONE}"
echo "============================================================"
