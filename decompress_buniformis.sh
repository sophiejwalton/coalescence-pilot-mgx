#!/bin/bash

#SBATCH --job-name=testgt-pro
#SBATCH --time=720
#SBATCH -p bigmem
#SBATCH -c 4
#SBATCH --mem=256GB

source /home/users/swalton/miniconda3/etc/profile.d/conda.sh

python3 workflow/scripts/track_snps_avg_v2_bootstrap_with_saving.py --outdir "workflow/report/track_snpsv2_ALL_bootstrapv3/" --indir "workflow/out/midas2_output/mergevfinal_101346/snps/" --species 101346