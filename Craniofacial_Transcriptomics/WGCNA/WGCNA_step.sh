#!/bin/bash

#SBATCH -J WGCNA_craniofacial
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=60G
#SBATCH --mail-type=end
#SBATCH --mail-user=yankee@uchc.edu


module load R/3.6.3

Rscript wgcna_step.R 
