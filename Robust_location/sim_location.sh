#!/bin/bash
#SBATCH --job-name=sim_robust_location
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_loc_%j.txt
#SBATCH --error=mc_error_loc_%j.txt
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --time=0-24:00:00
#SBATCH --mem=64G

# Load required cluster modules
module load gnu8 R gsl/2.6

# Prevent BLAS thread contention across the SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for Robust Location
# Default: 1000 MC reps, Welsh loss, n in {100, 500, 1000, 5000}, epsilon = 0.10, cv
Rscript simulation.R --reps=1000 --loss=Welsh --epsilon=0.10 --cv "$@"

# Generate final summary CSV and LaTeX table
Rscript read_sim_res.R

