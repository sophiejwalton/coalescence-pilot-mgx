#!/bin/bash

#SBATCH --job-name=testgt-pro
#SBATCH --time=72:00:00
#SBATCH -p dpetrov,hns
#SBATCH -c 6
#SBATCH --mem=48GB


source /home/users/swalton/miniconda3/etc/profile.d/conda.sh

conda activate /oak/stanford/groups/dpetrov/swalton/coalescence-pilot-mgx/.snakemake/conda/ea41405f604e8eeeead2c226ac256bd0_

midas2 merge_snps --samples_list workflow/out/list_of_samples.tsv --species_list 101346  --midasdb_name uhgg --robust_chunk --genome_depth 2 --site_depth 1  --site_prev 0.0 --snp_maf 0.01  --advanced --midasdb_dir /home/groups/bhgood/my_midasdb_uhgg --snp_type any --genome_coverage 0.5 --num_cores 6 workflow/out/midas2_output/mergevfinal_101346

