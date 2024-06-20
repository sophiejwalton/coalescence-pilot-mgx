import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
#from tqdm import tqdm
import argparse

from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
warnings.filterwarnings('ignore')

 
def get_main(species_dir, save_dir, species):
    info, depth, freq = load_and_sort_files(species_dir, species)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>5.]
    depth = depth[good_samples.index.values]

    freq = freq[good_samples.index.values]  
    depth_filtered= depth_filtering(depth)
    #depth_filtered.to_csv(f'{save_dir}/{species}_depth_filtered.csv.gz', compression='gzip')

    freq_filtered = freq_masked(freq, depth_filtered)
   # freq_filtered.to_csv(f'{save_dir}/{species}_freq_filtered.csv.gz',compression='gzip')
    #snps_info.to_csv(f'{save_dir}/{species}_snps_info.csv.gz',compression='gzip')
    num_int_sites, diversity = get_diversity_series(freq_filtered, thresh=.2)
    diversity_df = pd.DataFrame(diversity.reset_index()).rename(columns = {0:'diversity','index':'sample'})
    diversity_df.to_csv(f'{save_dir}/{species}_diversity_df.csv')
    med_nonzero_depth.to_csv(f'{save_dir}/{species}_med_depth_df.csv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species)
    


        
    

