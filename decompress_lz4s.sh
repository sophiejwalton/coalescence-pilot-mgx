#!/bin/bash

#SBATCH --job-name=testgt-pro
#SBATCH --time=720
#SBATCH -p bigmem
#SBATCH -c 4
#SBATCH --mem=256GB

source /home/users/swalton/miniconda3/etc/profile.d/conda.sh

lz4 -d workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_freqs.tsv.lz4 workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_freqs.tsv
#lz4 -d workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_depth.tsv.lz4 workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_depth.tsv 
#lz4 -d workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_info.tsv.lz4 workflow/out/midas2_output/mergevfinal3_101346/snps/101346/101346.snps_info.tsv
