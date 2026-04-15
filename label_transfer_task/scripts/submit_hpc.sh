#!/bin/bash
#SBATCH --job-name=label-transfer-benchmark
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=normal
#SBATCH --output=slurm-%j.log

# Label-Transfer Benchmark pipeline on HPC (SLURM)
#
# Usage:
#   sbatch scripts/submit_hpc.sh
#
# This submits the label-transfer benchmark pipeline with SLURM parallelization

# Load required modules (customize for your HPC)
module purge 2>/dev/null || true
module load GCC/12.3.0 2>/dev/null || true
module load R/4.3.2 2>/dev/null || true
module load GLPK/5.0 2>/dev/null || true
module load cairo/1.17.8 2>/dev/null || true
module load freetype/2.13.0 2>/dev/null || true
module load libwebp/1.3.1 2>/dev/null || true

cd "$(dirname "$(readlink -f "$0")")/.." || exit 1

echo "[$(date)] Starting label-transfer benchmark pipeline on HPC..."

# Run targets with SLURM future plan
Rscript - <<'EOF'
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2024-01-15"))
project_root <- normalizePath("..")

# Force project renv library on .libPaths()
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
if (dir.exists(renv_lib)) {
  .libPaths(unique(c(renv_lib, .libPaths())))
}

cat("[INFO] .libPaths():\n")
writeLines(.libPaths())

library(targets)
library(future.batchtools)

# Configure SLURM future plan
plan(
  future.batchtools::batchtools_slurm,
  workers = 12,  # Fewer workers than ISC (smaller jobs)
  resources = list(
    n_cores = 4,
    memory_gb = 32,
    walltime = "04:00:00",
    partition = "normal"
  )
)

# Run targets workflow
tar_make()

# Print summary
cat("\n")
tar_summary()
EOF

echo "[$(date)] Pipeline complete."
