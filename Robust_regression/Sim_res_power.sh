#!/bin/bash
#SBATCH --job-name=sim_tv_reg_power
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_power_%j.txt
#SBATCH --error=mc_error_power_%j.txt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --time=0-48:00:00
#SBATCH --mem=150G

# Load required cluster modules
module load gnu8 R gsl/2.6

# Prevent BLAS thread contention across the 80 SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for POWER (H1 Abrupt & H2 Gradual)
# Default: 1000 MC reps, 100 bootstrap draws, n in {500, 1000, 5000}
# Contamination: Evaluates all 7 regimes (Clean, AO 5/10%, IO 5/10%, RO 5/10%)
# Delta grid: c(0.00, 0.25, 0.50, 0.75, 1.00, 1.25) across Model I & Model II
# Outputs: 7 standalone figures (PNG/PDF), unified appendix PDF, and LaTeX tables
Rscript simulation_power.R "$@"


