#!/bin/bash
#SBATCH --job-name=sim_tv_reg_size
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_size_%j.txt
#SBATCH --error=mc_error_size_%j.txt
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --time=0-48:00:00
#SBATCH --mem=150G

# Load required cluster modules
module load gnu8 R gsl/2.6

# Prevent BLAS thread contention across the SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for empirical SIZE (Table 2)
# Default: 1000 MC reps, 100 bootstrap iterations, n in {500, 1000, 5000}
Rscript simulation_size.R "$@"

