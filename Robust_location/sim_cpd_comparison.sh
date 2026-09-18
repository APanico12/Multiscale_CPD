#!/bin/bash
#SBATCH --job-name=sim_cpd_comp
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_cpd_%j.txt
#SBATCH --error=mc_error_cpd_%j.txt
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=0-48:00:00
#SBATCH --mem=120G

# Load required cluster modules
module load gnu8 R gsl/2.6

# Note: Before submitting this job for the first time, ensure 'robcp' is installed
# in your personal R library. Execute once from the cluster login node:
# Rscript -e 'for (p in c("robcp", "doParallel", "foreach")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")'

# Prevent BLAS thread contention across the SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for CPD Comparison
# Default: 500 reps, n in {200, 500, 1000}, epsilon = 0.05, shift = 1.25
#Rscript simulation_cpd_comparison.R --reps=500 --n=200,500,1000 --epsilon=0.05 --shift=1.25 "$@"
# Default: 500 reps, n in {200, 500, 700}, epsilon = 0.05, shift = 0.5, var_scenarios = i,ii,iii
Rscript simulation_cpd_comparison.R --reps=500 --n=200,500,700 --epsilon=0.05 --shift=0.5 --var_scenarios=i,ii,iii "$@"

