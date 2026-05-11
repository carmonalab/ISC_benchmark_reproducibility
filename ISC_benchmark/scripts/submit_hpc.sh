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
  ensure_r_environment
  cd "$PROJECT_DIR/ISC_benchmark"

  Rscript --vanilla - <<'RS'
project_root <- normalizePath("..")
activate <- file.path(project_root, "renv", "activate.R")
if (file.exists(activate)) {
  source(activate)
}
if (requireNamespace("renv", quietly = TRUE)) {
  renv::load(project = project_root)
}

suppressPackageStartupMessages(library(yaml))

config <- yaml::read_yaml("config/dataset_idents.yaml")
available_ids <- names(config$idents)

processed_files <- list.files(
  file.path(project_root, "data", "processed"),
  pattern = "\\.rds$",
  full.names = FALSE
)
processed_ids <- unique(sub("_[^_]+\\.rds$", "", processed_files))

dataset_ids <- intersect(available_ids, processed_ids)

requested <- c(
  Sys.getenv("ISC_DATASET_IDS", unset = ""),
  Sys.getenv("ISC_DATASET_ID", unset = ""),
  Sys.getenv("ISC_TEST_DATASET", unset = "")
)
requested <- requested[nzchar(requested)]
if (length(requested) > 0) {
  requested <- unique(trimws(unlist(strsplit(paste(requested, collapse = ","), ",", fixed = TRUE))))
  missing_ids <- setdiff(requested, dataset_ids)
  if (length(missing_ids) > 0) {
    stop("Requested dataset(s) not available: ", paste(missing_ids, collapse = ", "))
  }
  dataset_ids <- requested
}

if (length(dataset_ids) == 0) {
  stop("No datasets found to submit. Check data/processed and config/dataset_idents.yaml.")
}

cat(paste(dataset_ids, collapse = "\n"))
RS
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
