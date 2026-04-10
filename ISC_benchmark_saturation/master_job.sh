#!/bin/bash -l

#SBATCH --account scarmona_til_omics
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user josep.garnicacaparros@unil.ch

#SBATCH --chdir /scratch/jgarnica/ISC_benchmark_reproducibility/ISC_benchmark_saturation
#SBATCH --job-name isc_saturation
#SBATCH --output /scratch/jgarnica/ISC_benchmark_reproducibility/reports/isc_saturation.out
#SBATCH --error /scratch/jgarnica/ISC_benchmark_reproducibility/reports/isc_saturation.err

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=2:00:00
#SBATCH --export=NONE

module load r-light/4.3.3

Rscript scripts/run_saturation_pipeline.R
