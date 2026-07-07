#!/bin/bash
# submit_full_pipeline_per_core.sh
# Run full CODEX core processing (segmentation-derived cell table -> tumor-nest boundary,
# myeloid infiltration, pEMT distance metrics) for every TMA core listed in cores.txt.
#
# Portability note: the original version of this script ran on an LSF cluster, submitting one
# independent `bsub` job per core inside an Apptainer/Singularity container (LSF resource flags
# -q/-n/-R/-M, environment-module loading, LSB_JOBTMP/LSB_JOBID-based scratch isolation). None of
# that is portable to an arbitrary system, so this version runs each core sequentially with a
# plain Rscript call instead. HPC users adapting this back to a scheduler should wrap the single
# Rscript call below in their own job-submission command (bsub/sbatch/etc.) with cluster-specific
# resource values of their own choosing.

set -euo pipefail

# -------------------- EDIT THESE PATHS --------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RUN_R="${REPO_ROOT}/analysis_scripts/codex_core_full_job.R"
TOOLS_R="${REPO_ROOT}/analysis_scripts/functions/CODEX_tools.R"
CELLS_CSV="${REPO_ROOT}/results/Generated_Data/CODEX/cell_type_and_state_assign_v2.csv"
MASK_DIR="${REPO_ROOT}/data/CODEX/segmentation_masks"
BOUND_DIR="${REPO_ROOT}/data/CODEX/tumor_boundary_masks"

OUT_DIR_RES="${REPO_ROOT}/results/Generated_Data/CODEX/processed_cores"
OUT_DIR_OUT="${REPO_ROOT}/results/Generated_Data/CODEX/per_core_results"
OUT_DIR_PLOT="${REPO_ROOT}/results/Paper_Figures/Revision_Figs/CODEX_TMA/Cores_plots"

CORES_LIST="${REPO_ROOT}/data/CODEX/cores.txt"

# Optional flip-map CSV with columns: core,flip (TRUE/FALSE). Leave empty to disable.
FLIP_MAP_CSV="${REPO_ROOT}/data/CODEX/flip_map.csv"

# Pixel size (µm per pixel)
PX2UM="0.5067"

# If TRUE: rebuild outputs even if Core_<ID>.rds already exists
FORCE_REBUILD="FALSE"

# Distance settings
SOURCE_LABEL="pEMT"
BOUNDARY="TRUE"
# Format: Name=Type1|Type2;Name2=Type3   (pipes separate types; semicolons separate groups)
TARGETS='Myeloid=Myeloid|Neutrophil|M2-like_Macrophage;Fibroblast=Fibroblast'

# Number of CPU cores used per-core for the internal distance computation (parallel::mclapply)
NCORES=4

# -------------------- END USER EDITS --------------------

mkdir -p "${OUT_DIR_RES}" "${OUT_DIR_OUT}" "${OUT_DIR_PLOT}"

FLIP_ARG=""
if [[ -n "${FLIP_MAP_CSV}" && -f "${FLIP_MAP_CSV}" ]]; then
  FLIP_ARG="flip_map_csv=\"${FLIP_MAP_CSV}\""
fi

# Read one core ID per line (ignore blank/comment lines, handle CRLF), process sequentially
while IFS= read -r CORE || [[ -n "${CORE:-}" ]]; do
  CORE="${CORE%$'\r'}"
  CORE="${CORE//[$'\t\r\n ']/}"
  [[ -z "${CORE}" || "${CORE:0:1}" == "#" ]] && continue

  RES_RDS="${OUT_DIR_RES}/Core_${CORE}.rds"
  OUT_RDS="${OUT_DIR_OUT}/Core_${CORE}.rds"

  if [[ -f "${RES_RDS}" && -f "${OUT_RDS}" && "${FORCE_REBUILD}" != "TRUE" ]]; then
    echo "Skipping ${CORE} (both outputs exist; FORCE_REBUILD=${FORCE_REBUILD})"
    continue
  fi

  echo "Processing core ${CORE}..."
  Rscript --vanilla "${RUN_R}" \
    core_id="${CORE}" \
    cells_csv="${CELLS_CSV}" \
    mask_dir="${MASK_DIR}" \
    bound_dir="${BOUND_DIR}" \
    tools_r="${TOOLS_R}" \
    out_dir_res="${OUT_DIR_RES}" \
    out_dir_out="${OUT_DIR_OUT}" \
    out_dir_plot="${OUT_DIR_PLOT}" \
    ${FLIP_ARG} \
    px2um="${PX2UM}" \
    source_label="${SOURCE_LABEL}" \
    boundary="${BOUNDARY}" \
    targets='${TARGETS}' \
    force_rebuild="${FORCE_REBUILD}" \
    ncores="${NCORES}"

  echo "Done: ${CORE}"
done < "${CORES_LIST}"
