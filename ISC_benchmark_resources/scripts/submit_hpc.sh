#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "${SCRIPT_DIR}/../.." && pwd )"
RESOURCES_DIR="${PROJECT_ROOT}/ISC_benchmark_resources"
LOG_DIR="${RESOURCES_DIR}/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

SLURM_PARTITION="${RESOURCE_SLURM_PARTITION:-public-cpu}"
SLURM_NODES="${RESOURCE_SLURM_NODES:-1}"
SLURM_CPUS="${RESOURCE_SLURM_CPUS:-4}"
SLURM_MEM="${RESOURCE_SLURM_MEM:-499G}"
SLURM_TIME="${RESOURCE_SLURM_TIME:-24:00:00}"

mkdir -p "${LOG_DIR}"

log_message() {
  local level="$1"
  shift
  printf '[%s] [%s] %s\n' "$(date +'%Y-%m-%d %H:%M:%S')" "${level}" "$*" >&2
}

trim_whitespace() {
  local value="$1"
  value="${value#${value%%[![:space:]]*}}"
  value="${value%${value##*[![:space:]]}}"
  printf '%s' "${value}"
}

read_dataset_ids() {
  cd "${RESOURCES_DIR}"

  local requested_raw=""
  local requested_value
  for requested_value in "${RESOURCE_DATASET_IDS:-}" "${RESOURCE_DATASET_ID:-}" "${RESOURCE_TEST_DATASET:-}"; do
    if [[ -n "${requested_value}" ]]; then
      requested_raw+="${requested_raw:+,}${requested_value}"
    fi
  done

  local -a processed_ids=()
  local -a selected_ids=()
  local -A family_lookup=()
  local -A processed_lookup=()
  local -A seen_requested=()

  while IFS= read -r dataset_id; do
    if [[ -n "${dataset_id}" ]]; then
      family_lookup["${dataset_id}"]=1
    fi
  done < <(
    awk '
      /^idents:/ { in_idents = 1; next }
      in_idents && /^[[:space:]]{2}[A-Za-z0-9_.-]+:[[:space:]]*$/ {
        gsub(/^[[:space:]]+/, "", $1)
        gsub(/:$/, "", $1)
        print $1
      }
    ' "${PROJECT_ROOT}/ISC_benchmark/config/dataset_idents.yaml"
  )

  while IFS= read -r processed_id; do
    [[ -z "${processed_id}" ]] && continue
    local_family="${processed_id%%_*}"

    if [[ -n "${family_lookup[${local_family}]:-}" ]]; then
      processed_lookup["${processed_id}"]=1
      processed_ids+=("${processed_id}")
    fi
  done < <(
    find "${PROJECT_ROOT}/data/processed" -maxdepth 1 -name '*.rds' -exec basename {} .rds \; | sort -u
  )

  if [[ -n "${requested_raw}" ]]; then
    IFS=',' read -r -a requested_ids <<< "${requested_raw}"
    for dataset_id in "${requested_ids[@]}"; do
      dataset_id="$(trim_whitespace "${dataset_id}")"
      [[ -z "${dataset_id}" ]] && continue
      if [[ -z "${processed_lookup[${dataset_id}]:-}" ]]; then
        log_message "ERROR" "Requested dataset is not available in data/processed or not configured: ${dataset_id}"
        exit 1
      fi
      if [[ -z "${seen_requested[${dataset_id}]:-}" ]]; then
        selected_ids+=("${dataset_id}")
        seen_requested["${dataset_id}"]=1
      fi
    done
  else
    selected_ids=("${processed_ids[@]}")
  fi

  if [[ ${#selected_ids[@]} -eq 0 ]]; then
    log_message "ERROR" "No datasets found to submit."
    exit 1
  fi

  printf '%s\n' "${selected_ids[@]}"
}

mapfile -t DATASET_IDS < <(read_dataset_ids)

for dataset_id in "${DATASET_IDS[@]}"; do
  [[ -z "${dataset_id}" ]] && continue

  log_message "INFO" "Submitting resource benchmark for ${dataset_id}"

  sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=isc_resources_${dataset_id}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME}
#SBATCH --nodes=${SLURM_NODES}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --output=${LOG_DIR}/${dataset_id}_${TIMESTAMP}_%j.out
#SBATCH --error=${LOG_DIR}/${dataset_id}_${TIMESTAMP}_%j.err

export RESOURCE_DATASET_ID="${dataset_id}"
export RESOURCE_SLURM_CPUS="${SLURM_CPUS}"

bash "${RESOURCES_DIR}/scripts/master_job.sh"
EOF
done