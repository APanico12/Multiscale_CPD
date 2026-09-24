#!/bin/bash
#SBATCH --job-name=cpd_sigma_ii
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_cpd_sigma_ii_%j.txt
#SBATCH --error=mc_error_cpd_sigma_ii_%j.txt
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=0-48:00:00
#SBATCH --mem=120G

# ==============================================================================
# AVAILABLE CHANGE-POINT DETECTION MODELS / TESTS IN COMPARISON:
# ------------------------------------------------------------------------------
#   W   : Our Proposed Linearized CUSUM test with Welsh score (Our_Welsh)
#   L   : Two-sample Hodges-Lehmann test (Hodges_Lehmann via robcp::hl_test) [also 'HL']
#   H   : Huberized CUSUM test (Huber_CUSUM via robcp::huber_cusum)
#   WMW : Wilcoxon-Mann-Whitney rank test (Wilcoxon via robcp::wmw_test)
#   S   : Schmidt (2021) Gini test for heteroscedastic time series (Schmidt_Gini)
# ==============================================================================

# --- USER-CONFIGURABLE PARAMETERS ---

# 1. Models to include in the simulation:
#    Select any combination: e.g. "L,W,S" or "W,L,H,WMW,S" or "W,S"
MODELS="W,H,WMW,S"

# 2. Bandwidth rate exponent k (e.g. 0.45 -> window size k_n = floor(n^0.45)):
#    Used when USE_CV="FALSE".
BANDWIDTH="0.45"

# 3. Cross-Validation (CV) for bandwidth selection (Welsh test):
#    Set USE_CV="TRUE" to select k_n per replication via forward predictive CV,
#    or USE_CV="FALSE" to use the fixed rate specified in BANDWIDTH above.
USE_CV="FALSE"
CV_GRID="0.45,0.65"

# 4. Monte Carlo replications and design grid:
REPS=500
N_VALS="500,1000,5000"
EPSILON="0.05"
SHIFT="0.5"
VAR_SCENARIOS="ii"

# 5. Output CSV file names:
OUTPUT_RAW="sim_results_cpd_comparison_sigma_ii.csv"
OUTPUT_SUMMARY="sim_summary_cpd_comparison_sigma_ii.csv"

# Load required cluster modules
module load gnu8 R gsl/2.6

# Prevent BLAS thread contention across the SLURM tasks
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute parallel Monte Carlo simulation for CPD Comparison (Variance Scenario ii)
Rscript simulation_cpd_comparison.R \
  --reps=${REPS} \
  --n=${N_VALS} \
  --epsilon=${EPSILON} \
  --shift=${SHIFT} \
  --var_scenarios=${VAR_SCENARIOS} \
  --k=${BANDWIDTH} \
  --models="${MODELS}" \
  --use_cv=${USE_CV} \
  --cv_grid=${CV_GRID} \
  --output_file="${OUTPUT_RAW}" \
  --output_summary="${OUTPUT_SUMMARY}" \
  "$@"

