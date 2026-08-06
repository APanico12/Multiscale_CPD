#!/bin/bash
#SBATCH --job-name=sim_alt_clayton_linear_d_2_tau_m2
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_%j.txt
#SBATCH --error=mc_error_%j.txt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --time=0-48:00:00
#SBATCH --mem=350G
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL

# Load the required modules
module load gnu8 R gsl/2.6

# Run the R script
Rscript sim_alt_clayton_linear_d_2_tau_m2.R
