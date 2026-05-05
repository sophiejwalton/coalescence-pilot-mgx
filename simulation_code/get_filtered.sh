#!/bin/bash

#SBATCH --job-name=get_filtered
#SBATCH --time=4:00:00
#SBATCH -p bigmem
#SBATCH --mem=256GB
#SBATCH -c 4


source /home/users/swalton/miniconda3/etc/profile.d/conda.sh

python3 workflow/scripts/track_snps_avg_v2_bootstrap_iss_v2.py --outdir workflow/report/track_snps/ --indir workflow/out/midas2_output20/merge_100196/snps --species 100196
