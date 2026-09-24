#!/bin/bash
#SBATCH --job-name=sim_tv_reg_power_LMUH
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_power_LMUH_%j.txt
#SBATCH --error=mc_error_power_LMUH_%j.txt
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --time=0-48:00:00
#SBATCH --mem=150G

# --- USER-CONFIGURABLE PARAMETERS ---
# 1. Bandwidth rate exponent k (e.g. 0.45 -> window size k_n = floor(n^0.45)):
#    Used when USE_CV="FALSE".
BANDWIDTH="0.45"

# 2. Cross-Validation (CV) for bandwidth selection:
#    Set USE_CV="TRUE" to select k_n per replication via forward predictive CV,
#    or USE_CV="FALSE" to use the fixed rate specified in BANDWIDTH above.
USE_CV="FALSE"
CV_GRID="0.45,0.65"

# 3. Monte Carlo replications and design grid:
REPS=1000
B_BOOT=100
N_VALS="500,1000,5000"
DELTA="0.0,0.25,0.50,0.75,1.00,1.25"

# Load required cluster modules
module load gnu8 R gsl/2.6

# Prevent BLAS thread contention across the SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for POWER (H1 Abrupt & H2 Gradual)
# Model: Model II (LMUH: Unconditional Heteroscedasticity)
Rscript simulation_power.R \
  --LMUH \
  --reps=${REPS} \
  --boot=${B_BOOT} \
  --n=${N_VALS} \
  --delta="${DELTA}" \
  --k=${BANDWIDTH} \
  --use_cv=${USE_CV} \
  --cv_grid=${CV_GRID} \
  "$@"


