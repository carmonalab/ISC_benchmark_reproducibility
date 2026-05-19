#!/bin/bash
set -euo pipefail

# ISC Benchmark pipeline on HPC (SLURM)
#
# Usage:
#   bash scripts/submit_hpc.sh
#
# Dataset selection:
#   ISC_DATASET_ID=JoaI_CRC-SG1_Normal bash scripts/submit_hpc.sh
#   ISC_TEST_DATASET=JoaI_CRC-SG1_Normal bash scripts/submit_hpc.sh
#   ISC_DATASET_IDS="JoaI_CRC-SG1_Normal,Stephenson_Cambridge_Covid" bash scripts/submit_hpc.sh
#
# The wrapper submits one SLURM job per dataset. Each job runs the targets
# pipeline with that dataset filtered at the graph-construction stage.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
REPO_DIR="$( git -C "$SCRIPT_DIR" rev-parse --show-toplevel )"
LOG_DIR="$PROJECT_DIR/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

SLURM_PARTITION="${SLURM_PARTITION:-public-bigmem}"
SLURM_NODES="${SLURM_NODES:-1}"
SLURM_CPUS="${SLURM_CPUS:-4}"
SLURM_MEM="${SLURM_MEM:-950G}"
SLURM_TIME="${SLURM_TIME:-48:00:00}"

DATASET_TASKS=(missclassify SplitCelltype Nct cellular_complexity Nsamples NCell)
CROSS_DATASET_TASKS=(batch_effects biological_perturbations)

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
  cd "$PROJECT_DIR"

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

    if [[ -n "${family_lookup[$local_family]:-}" ]]; then
      processed_lookup["$processed_id"]=1
      processed_ids+=("$processed_id")
    fi
  done < <(
    find "$REPO_DIR/data/processed" -maxdepth 1 -name '*.rds' -printf '%f\n' | sed -E 's/\.rds$//' | sort -u
  )

  if [[ -n "$requested_raw" ]]; then
    IFS=',' read -r -a requested_ids <<< "$requested_raw"
    for dataset_id in "${requested_ids[@]}"; do
      dataset_id="$(trim_whitespace "$dataset_id")"
      [[ -z "$dataset_id" ]] && continue
      if [[ -z "${processed_lookup[$dataset_id]:-}" ]]; then
        log_message "ERROR" "Requested dataset stem is not available in repo data/processed or not configured: $dataset_id"
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
    log_message "ERROR" "No datasets found to submit. Check repo data/processed and config/dataset_idents.yaml."
    exit 1
  fi

  printf '%s\n' "${selected_ids[@]}"
}

read_dataset_families() {
  cd "$PROJECT_DIR"

  local requested_raw=""
  local requested_value
  for requested_value in "${ISC_DATASET_FAMILIES:-}" "${ISC_DATASET_FAMILY:-}"; do
    if [[ -n "$requested_value" ]]; then
      requested_raw+="${requested_raw:+,}${requested_value}"
    fi
  done
  local -a families=()
  local -A family_lookup=()
  local -A requested_lookup=()
  local has_family_filter=false

  if [[ -n "$requested_raw" ]]; then
    has_family_filter=true
    IFS=',' read -r -a requested_families <<< "$requested_raw"
    for family_name in "${requested_families[@]}"; do
      family_name="$(trim_whitespace "$family_name")"
      [[ -z "$family_name" ]] && continue
      requested_lookup["$family_name"]=1
    done
  fi

  while IFS= read -r dataset_stem; do
    [[ -z "$dataset_stem" ]] && continue
    family_name="${dataset_stem%%_*}"

    if [[ "$has_family_filter" == true && -z "${requested_lookup[$family_name]:-}" ]]; then
      continue
    fi

    if [[ -z "${family_lookup[$family_name]:-}" ]]; then
      families+=("$family_name")
      family_lookup["$family_name"]=1
    fi
  done < <(
    find "$REPO_DIR/data/processed" -maxdepth 1 -name '*.rds' -printf '%f\n' | sed -E 's/\.rds$//' | sort -u
  )

  if [[ ${#families[@]} -eq 0 ]]; then
    log_message "ERROR" "No dataset families found in repo data/processed"
    exit 1
  fi

  printf '%s\n' "${families[@]}"
}

run_dataset_locally() {
  local dataset_id="$1"
  local tasks_csv="$2"
  log_message "INFO" "Running dataset locally: $dataset_id"

  ensure_r_environment
  cd "$PROJECT_DIR"
  ISC_DATASET_ID="$dataset_id" ISC_TASKS="$tasks_csv" Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
project_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(project_lib)) {
  .libPaths(unique(c(project_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
cat("[INFO] ISC_DATASET_ID = ", Sys.getenv("ISC_DATASET_ID"), "\n", sep = "")

library(targets)

selector <- Sys.getenv("ISC_DATASET_ID", unset = "")
if (!nzchar(selector)) {
  selector <- Sys.getenv("ISC_DATASET_FAMILIES", unset = "")
}
if (!nzchar(selector)) {
  selector <- "default"
}
safe_selector <- gsub("[^A-Za-z0-9._-]", "_", selector)
store_dir <- file.path("_targets", paste0("store_", safe_selector))
dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] targets store = ", store_dir, "\n", sep = "")

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL, store = store_dir)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS
}

run_family_locally() {
  local family_id="$1"
  local tasks_csv="$2"
  log_message "INFO" "Running family locally: $family_id"

  ensure_r_environment
  cd "$PROJECT_DIR"
  ISC_DATASET_FAMILIES="$family_id" ISC_TASKS="$tasks_csv" Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
project_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(project_lib)) {
  .libPaths(unique(c(project_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
cat("[INFO] ISC_DATASET_FAMILIES = ", Sys.getenv("ISC_DATASET_FAMILIES"), "\n", sep = "")
cat("[INFO] ISC_TASKS = ", Sys.getenv("ISC_TASKS"), "\n", sep = "")

library(targets)

selector <- Sys.getenv("ISC_DATASET_ID", unset = "")
if (!nzchar(selector)) {
  selector <- Sys.getenv("ISC_DATASET_FAMILIES", unset = "")
}
if (!nzchar(selector)) {
  selector <- "default"
}
safe_selector <- gsub("[^A-Za-z0-9._-]", "_", selector)
store_dir <- file.path("_targets", paste0("store_", safe_selector))
dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] targets store = ", store_dir, "\n", sep = "")

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL, store = store_dir)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS
}

submit_targets_job() {
  local job_type="$1"
  local job_label="$2"
  local tasks_csv="$3"
  local selector_name="$4"
  local selector_value="$5"
  local safe_label
  safe_label="$(printf '%s' "$job_label" | tr -c 'A-Za-z0-9._-' '_')"

  log_message "INFO" "Submitting ${job_type} job: ${job_label}"

  local -a env_exports=(
    "PROJECT_DIR=$PROJECT_DIR"
    "ISC_TASKS=$tasks_csv"
    "ISC_DATASET_IDS="
    "ISC_DATASET_ID="
    "ISC_TEST_DATASET="
    "ISC_DATASET_FAMILIES="
    "ISC_DATASET_FAMILY="
    "${selector_name}=$selector_value"
  )

  env "${env_exports[@]}" sbatch --parsable \
    --export=ALL \
    --job-name="isc-${safe_label}" \
    --partition="$SLURM_PARTITION" \
    --nodes="$SLURM_NODES" \
    --cpus-per-task="$SLURM_CPUS" \
    --mem="$SLURM_MEM" \
    --time="$SLURM_TIME" \
    --output="$LOG_DIR/slurm-${safe_label}-%j.out" \
    --error="$LOG_DIR/slurm-${safe_label}-%j.err" <<'SLURM_SCRIPT'
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

cd "$PROJECT_DIR"

echo "========== Job Started =========="
echo "Time: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "ISC_DATASET_ID: ${ISC_DATASET_ID:-unset}"
echo "ISC_DATASET_FAMILIES: ${ISC_DATASET_FAMILIES:-unset}"
echo "ISC_TASKS: ${ISC_TASKS:-unset}"
echo "================================="

Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
project_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(project_lib)) {
  .libPaths(unique(c(project_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())
cat("[INFO] ISC_DATASET_ID = ", Sys.getenv("ISC_DATASET_ID"), "\n", sep = "")
cat("[INFO] ISC_DATASET_FAMILIES = ", Sys.getenv("ISC_DATASET_FAMILIES"), "\n", sep = "")
cat("[INFO] ISC_TASKS = ", Sys.getenv("ISC_TASKS"), "\n", sep = "")

library(targets)

selector <- Sys.getenv("ISC_DATASET_ID", unset = "")
if (!nzchar(selector)) {
  selector <- Sys.getenv("ISC_DATASET_FAMILIES", unset = "")
}
if (!nzchar(selector)) {
  selector <- "default"
}
safe_selector <- gsub("[^A-Za-z0-9._-]", "_", selector)
store_dir <- file.path("_targets", paste0("store_", safe_selector))
dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] targets store = ", store_dir, "\n", sep = "")

message("Starting targets::tar_make(callr_function = NULL) at ", Sys.time())
targets::tar_make(callr_function = NULL, store = store_dir)
message("Finished targets::tar_make(callr_function = NULL) at ", Sys.time())
RS

echo "========== Job Finished =========="
echo "Time: $(date)"
echo "================================="
SLURM_SCRIPT
}

submit_dataset_job() {
  local dataset_id="$1"
  local tasks_csv="$2"
  submit_targets_job "dataset" "$dataset_id" "$tasks_csv" "ISC_DATASET_ID" "$dataset_id"
}

submit_family_job() {
  local family_id="$1"
  local tasks_csv="$2"
  submit_targets_job "family" "$family_id" "$tasks_csv" "ISC_DATASET_FAMILIES" "$family_id"
}

log_message "INFO" "========== ISC Benchmark HPC Submission =========="
log_message "INFO" "Timestamp: $TIMESTAMP"
log_message "INFO" "Project directory: $PROJECT_DIR"
log_message "INFO" "Log directory: $LOG_DIR"
log_message "INFO" "SLURM configuration: partition=$SLURM_PARTITION nodes=$SLURM_NODES cpus=$SLURM_CPUS mem=$SLURM_MEM time=$SLURM_TIME"

mapfile -t DATASET_IDS < <(read_dataset_ids)
DATASET_TASKS_CSV="$(IFS=,; echo "${DATASET_TASKS[*]}")"
CROSS_TASKS_CSV="$(IFS=,; echo "${CROSS_DATASET_TASKS[*]}")"
STEM_FILTER_SET="${ISC_DATASET_IDS:-}${ISC_DATASET_ID:-}${ISC_TEST_DATASET:-}"
FAMILY_FILTER_SET="${ISC_DATASET_FAMILIES:-}${ISC_DATASET_FAMILY:-}"
RUN_STEM_JOBS=true
RUN_FAMILY_JOBS=true

if [[ -n "$STEM_FILTER_SET" ]]; then
  RUN_FAMILY_JOBS=false
elif [[ -n "$FAMILY_FILTER_SET" ]]; then
  RUN_STEM_JOBS=false
fi

log_message "INFO" "Dataset selection: ${#DATASET_IDS[@]} stem job(s)"

if command -v sbatch >/dev/null 2>&1; then
  JOB_IDS=()

  if [[ "$RUN_STEM_JOBS" == true ]]; then
    for dataset_id in "${DATASET_IDS[@]}"; do
      job_id="$(submit_dataset_job "$dataset_id" "$DATASET_TASKS_CSV")"
      JOB_IDS+=("$job_id")
      log_message "INFO" "Submitted stem $dataset_id as SLURM job $job_id"
    done
  fi

  if [[ "$RUN_FAMILY_JOBS" == true ]]; then
    mapfile -t DATASET_FAMILIES < <(read_dataset_families)
    log_message "INFO" "Family selection: ${#DATASET_FAMILIES[@]} cross-dataset job(s)"
    for family_id in "${DATASET_FAMILIES[@]}"; do
      job_id="$(submit_family_job "$family_id" "$CROSS_TASKS_CSV")"
      JOB_IDS+=("$job_id")
      log_message "INFO" "Submitted family $family_id as SLURM job $job_id"
    done
  fi

  log_message "INFO" "Submission complete: ${JOB_IDS[*]}"
else
  if [[ "$RUN_STEM_JOBS" == true ]]; then
    log_message "WARNING" "sbatch not found. Running selected stem job(s) locally instead."
    for dataset_id in "${DATASET_IDS[@]}"; do
      run_dataset_locally "$dataset_id" "$DATASET_TASKS_CSV"
    done
  fi
  if [[ "$RUN_FAMILY_JOBS" == true ]]; then
    mapfile -t DATASET_FAMILIES < <(read_dataset_families)
    log_message "WARNING" "sbatch not found. Running family job(s) locally instead."
    for family_id in "${DATASET_FAMILIES[@]}"; do
      run_family_locally "$family_id" "$CROSS_TASKS_CSV"
    done
  fi
fi
