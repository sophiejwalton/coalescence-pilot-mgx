#!/bin/bash

#SBATCH --job-name=testgt-pro
#SBATCH --time=48:00:00
#SBATCH -p normal,hns
#SBATCH -c 6
#SBATCH --mem=64GB

python3 parse_data_SJW.py