#!/bin/bash
#SBATCH --job-name=sim_null_clayton_d_5_tau1_5_sine
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_%j.txt
#SBATCH --error=mc_error_%j.txt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60
#SBATCH --time=0-24:00:00
#SBATCH --mem=125G

# Load the required modules for the R environment
module load gnu8 R gsl/2.6

# Run the R script.
# The script is designed to automatically detect the 60 tasks allocated by SLURM.
echo "Starting simulation for sim_null_clayton_d_5_tau1_5_sine..."
Rscript sim_null_clayton_d_5_tau1_5_sine.R
echo "Simulation finished."
