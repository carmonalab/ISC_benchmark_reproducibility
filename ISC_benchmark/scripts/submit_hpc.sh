#!/bin/bash
set -euo pipefail

# ISC Benchmark pipeline on HPC (SLURM)
#
# Usage:
#   bash scripts/submit_hpc.sh
#
# Dataset selection:
#   ISC_DATASET_ID=JoaI bash scripts/submit_hpc.sh
#   ISC_TEST_DATASET=JoaI bash scripts/submit_hpc.sh
#   ISC_DATASET_IDS="JoaI,Stephenson" bash scripts/submit_hpc.sh
#
# The wrapper submits one SLURM job per dataset. Each job runs the targets
# pipeline with that dataset filtered at the graph-construction stage.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
LOG_DIR="$SCRIPT_DIR/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

SLURM_PARTITION="${SLURM_PARTITION:-normal}"
SLURM_NODES="${SLURM_NODES:-1}"
SLURM_CPUS="${SLURM_CPUS:-8}"
SLURM_MEM="${SLURM_MEM:-64G}"
SLURM_TIME="${SLURM_TIME:-24:00:00}"

mkdir -p "$LOG_DIR"

log_message() {
  local level="$1"
  shift
  printf '[%s] [%s] %s\n' "$(date +'%Y-%m-%d %H:%M:%S')" "$level" "$*" >&2
}

trim_whitespace() {
  local value="$1"
  value="${value#${value%%[![:space:]]*}}"
  value="${value%${value##*[![:space:]]}}"
  printf '%s' "$value"
}

ensure_r_environment() {
  if command -v Rscript >/dev/null 2>&1; then
    return 0
  fi

  if command -v module >/dev/null 2>&1; then
    module purge 2>/dev/null || true
    module load GCC/12.3.0 2>/dev/null || true
    module load R/4.3.2 2>/dev/null || true
    module load GLPK/5.0 2>/dev/null || true
    module load cairo/1.17.8 2>/dev/null || true
    module load freetype/2.13.0 2>/dev/null || true
    module load libwebp/1.3.1 2>/dev/null || true
  fi

  if ! command -v Rscript >/dev/null 2>&1; then
    log_message "ERROR" "Rscript is not available. Load the R module or configure PATH before running this script."
    exit 1
  fi
}

read_dataset_ids() {
  cd "$PROJECT_DIR/ISC_benchmark"

  local requested_raw=""
  local requested_value
  for requested_value in "${ISC_DATASET_IDS:-}" "${ISC_DATASET_ID:-}" "${ISC_TEST_DATASET:-}"; do
    if [[ -n "$requested_value" ]]; then
      requested_raw+="${requested_raw:+,}${requested_value}"
    fi
  done
  local -a processed_ids=()
  local -a selected_ids=()
  local -A family_lookup=()
  local -A processed_lookup=()
  local -A seen_requested=()

  while IFS= read -r dataset_id; do
    if [[ -n "$dataset_id" ]]; then
      family_lookup["$dataset_id"]=1
    fi
  done < <(
    awk '
      /^idents:/ { in_idents = 1; next }
      in_idents && /^[[:space:]]{2}[A-Za-z0-9_.-]+:[[:space:]]*$/ {
        gsub(/^[[:space:]]+/, "", $1)
        gsub(/:$/, "", $1)
        print $1
      }
    ' config/dataset_idents.yaml
  )

  while IFS= read -r processed_id; do
    [[ -z "$processed_id" ]] && continue

    local_family="${processed_id%%_*}"
    if [[ "$local_family" == "StephensonE" ]]; then
      local_family="Stephenson"
    fi

    if [[ -n "${family_lookup[$local_family]:-}" ]]; then
      processed_lookup["$processed_id"]=1
      processed_ids+=("$processed_id")
    fi
  done < <(
    find "$PROJECT_DIR/data/processed" -maxdepth 1 -name '*.rds' -printf '%f\n' | sed -E 's/\.rds$//' | sort -u
  )

  if [[ -n "$requested_raw" ]]; then
    IFS=',' read -r -a requested_ids <<< "$requested_raw"
    for dataset_id in "${requested_ids[@]}"; do
      dataset_id="$(trim_whitespace "$dataset_id")"
      [[ -z "$dataset_id" ]] && continue
      if [[ -z "${processed_lookup[$dataset_id]:-}" ]]; then
        log_message "ERROR" "Requested dataset stem is not available in data/processed or not configured: $dataset_id"
        exit 1
      fi
      if [[ -z "${seen_requested[$dataset_id]:-}" ]]; then
        selected_ids+=("$dataset_id")
        seen_requested["$dataset_id"]=1
      fi
    done
  else
    selected_ids=("${processed_ids[@]}")
  fi

  if [[ ${#selected_ids[@]} -eq 0 ]]; then
    log_message "ERROR" "No datasets found to submit. Check data/processed and config/dataset_idents.yaml."
    exit 1
  fi

  printf '%s\n' "${selected_ids[@]}"
}

run_dataset_locally() {
  local dataset_id="$1"
  log_message "INFO" "Running dataset locally: $dataset_id"

  ensure_r_environment
  cd "$PROJECT_DIR/ISC_benchmark"
  ISC_DATASET_ID="$dataset_id" Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) {
  source(activate)
}
if (requireNamespace("renv", quietly = TRUE)) {
  renv::load(project = project_root)
}

r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
cat("[INFO] ISC_DATASET_ID = ", Sys.getenv("ISC_DATASET_ID"), "\n", sep = "")

library(targets)

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS
}

submit_dataset_job() {
  local dataset_id="$1"
  local safe_dataset_id
  safe_dataset_id="$(printf '%s' "$dataset_id" | tr -c 'A-Za-z0-9._-' '_')"

  log_message "INFO" "Submitting dataset job: $dataset_id"

  sbatch --parsable \
    --job-name="isc-${safe_dataset_id}" \
    --partition="$SLURM_PARTITION" \
    --nodes="$SLURM_NODES" \
    --cpus-per-task="$SLURM_CPUS" \
    --mem="$SLURM_MEM" \
    --time="$SLURM_TIME" \
    --output="$LOG_DIR/slurm-${safe_dataset_id}-%j.out" \
    --error="$LOG_DIR/slurm-${safe_dataset_id}-%j.err" \
    --export=ALL,PROJECT_DIR="$PROJECT_DIR",ISC_DATASET_ID="$dataset_id" <<'SLURM_SCRIPT'
#!/bin/bash
set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  module load GCC/12.3.0 2>/dev/null || true
  module load R/4.3.2 2>/dev/null || true
  module load GLPK/5.0 2>/dev/null || true
  module load cairo/1.17.8 2>/dev/null || true
  module load freetype/2.13.0 2>/dev/null || true
  module load libwebp/1.3.1 2>/dev/null || true
fi

cd "$PROJECT_DIR/ISC_benchmark"

echo "========== Job Started =========="
echo "Time: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "ISC_DATASET_ID: ${ISC_DATASET_ID:-unset}"
echo "================================="

Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) {
  source(activate)
}
if (requireNamespace("renv", quietly = TRUE)) {
  renv::load(project = project_root)
}

r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
cat("[INFO] ISC_DATASET_ID = ", Sys.getenv("ISC_DATASET_ID"), "\n", sep = "")

library(targets)

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS

echo "========== Job Finished =========="
echo "Time: $(date)"
echo "================================="
SLURM_SCRIPT
}

log_message "INFO" "========== ISC Benchmark HPC Submission =========="
log_message "INFO" "Timestamp: $TIMESTAMP"
log_message "INFO" "Project directory: $PROJECT_DIR"
log_message "INFO" "Log directory: $LOG_DIR"
log_message "INFO" "SLURM configuration: partition=$SLURM_PARTITION nodes=$SLURM_NODES cpus=$SLURM_CPUS mem=$SLURM_MEM time=$SLURM_TIME"

mapfile -t DATASET_IDS < <(read_dataset_ids)
log_message "INFO" "Dataset selection: ${#DATASET_IDS[@]} dataset(s)"

if command -v sbatch >/dev/null 2>&1; then
  JOB_IDS=()
  for dataset_id in "${DATASET_IDS[@]}"; do
    job_id="$(submit_dataset_job "$dataset_id")"
    JOB_IDS+=("$job_id")
    log_message "INFO" "Submitted $dataset_id as SLURM job $job_id"
  done

  log_message "INFO" "Submission complete: ${JOB_IDS[*]}"
else
  log_message "WARNING" "sbatch not found. Running selected dataset(s) locally instead."
  for dataset_id in "${DATASET_IDS[@]}"; do
    run_dataset_locally "$dataset_id"
  done
fi
