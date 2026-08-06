#!/bin/bash
#SBATCH --job-name=sim_null_clayton_d_2_tau1_5_const
#SBATCH --account=T_STAGE_LUIGI_GROSSI
#SBATCH --mail-user=antonio.panico@unipr.it
#SBATCH --mail-type=END,FAIL 
#SBATCH --output=mc_output_%j.txt
#SBATCH --error=mc_error_%j.txt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-24:00:00
#SBATCH --mem=125G



# Load the required modules for the R environment
module load gnu8 R gsl/2.6

# Run the R script.
echo "Starting simulation for sim_null_clayton_d_2_tau1_5_const..."
Rscript sim_null_clayton_d_2_tau1_5_const.R
echo "Simulation finished."