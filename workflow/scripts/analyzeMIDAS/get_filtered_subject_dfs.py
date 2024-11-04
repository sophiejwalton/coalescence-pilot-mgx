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

 
def get_main(species_dir, save_dir, species, genome_length):
    info, depth, freq = load_and_sort_files(species_dir)
    info = info.set_index('site_id')
    freq = freq.set_index('site_id')
    depth = depth.set_index('site_id')
    sites_considered = int(len(depth))
    depth_filtered = depth_filtering(depth)
    samples = freq.columns.values
    # get list of subjects
    subjects = []
    for sample in samples:
        subject = sample.split('-')[2]
        if subject not in subjects:
            subjects.append(subject)
    all_diversity = []
    all_depths  = []
    for i, subject in enumerate(subjects):

        subject_freq, subject_depth = get_subject_dfs(freq, depth_filtered, subject, samples)
        diversity_series, _= get_diversity_series(subject_freq, subject_depth, genome_length, sites_considered)
        median_depth_series = subject_depth.copy().replace(0, np.nan).median(skipna = True)
        good_depth_samples = median_depth_series[median_depth_series > 10].index.values
        good_depth_samples_good = list(good_depth_samples)

        if len(good_depth_samples_good) == 0:
            continue 
        good_freq, good_samples = cleanup_and_polarize(subject_freq[good_depth_samples_good], median_depth_series, diversity_series, subject)
        
        good_depth = subject_depth[good_samples]
        good_freq = good_freq[good_samples]
        good_depth_na = good_depth*good_depth.isna() + 1. 
        good_freq = good_freq*good_depth_na
        
        good_depth, good_freq = filter_sites_across_samples(good_depth, good_freq)

        good_depth.to_csv(f'{save_dir}/{subject}_good_depth.csv.gz',compression='gzip')
        good_freq.to_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip')

        
        diversity_series, _ = get_diversity_series(good_freq, good_depth, genome_length, sites_considered)
        
        median_depth_series = median_depth_series[good_freq.columns.values]
        median_depth_series = pd.DataFrame(median_depth_series).reset_index()
        
        
        median_depth_series['timepoint'] = median_depth_series['index'].str.split("-").str[-1]
        median_depth_series['timepoint'] = pd.to_numeric(median_depth_series['timepoint'])
        diversity_series = pd.DataFrame(diversity_series).reset_index()
        diversity_series['timepoint'] = diversity_series['index'].str.split("-").str[-1]
        diversity_series['timepoint'] = pd.to_numeric(diversity_series['timepoint'])
        diversity_series['Subject'] = subject
        diversity_series['Species'] = species
        median_depth_series['Subject'] = subject
        median_depth_series['Species'] = species
        all_diversity.append(diversity_series)
        all_depths.append(median_depth_series)
    all_diversity = pd.concat(all_diversity)
    all_depths = pd.concat(all_depths)
    all_diversity.to_csv(f'{save_dir}/diversity_df.csv')
    all_depths.to_csv(f'{save_dir}/median_depth_df.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter sites across samples')

    # add arguments
    parser.add_argument('out_dir', action='store',
                    help='Outdir prefix where files are stored')
    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')
    parser.add_argument('species', action = 'store', 
                       help = 'species to perform analysis on')
    args = parser.parse_args()
    species_info = pd.read_csv('workflow/reports/species_info_full.csv')
    genome_length = species_info.loc[species_info['species_id'] == args.species, 'length'].values[0]
    species_dir = f'{args.out_dir}/{args.species}'
    save_dir = f'{args.save_dir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)

    get_main(species_dir, save_dir, args.species, genome_length)
    


        
    

