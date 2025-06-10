import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
import itertools as it
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
warnings.filterwarnings('ignore')


def get_fixed_diffs(freq_filtered,sample):
    freq_filtered_sample = freq_filtered[[sample]].copy()

  #  good_sites_lower = freq_filtered_sample.loc[freq_filtered_sample[sample] < .2].index.values
    fixed_diffs_lower_snps = freq_filtered.loc[freq_filtered[sample] < .2,:] >.8 
    fixed_diffs_lower_snps = fixed_diffs_lower_snps.sum(axis=0)

   # good_sites_upper = freq_filtered_sample.loc[freq_filtered_sample[sample] > .8].index.values
    fixed_diffs_upper_snps = freq_filtered.loc[freq_filtered[sample] >.8,:] <.2
    fixed_diffs_upper_snps = fixed_diffs_upper_snps.sum(axis=0)

    num_sites_non_int = (freq_filtered >.8).sum(axis=0) + (freq_filtered <.2).sum(axis=0) 
    fixed_diffs = fixed_diffs_upper_snps+fixed_diffs_lower_snps

    fixed_diffs = fixed_diffs.rename('fixed_diffs')
    num_sites_non_int = num_sites_non_int.rename('comparisons')
 
    return pd.concat([fixed_diffs,num_sites_non_int],axis=1)
 
def get_main(species_dir, save_dir, species,):
    metadata = pd.read_csv('config/e003_metadata_cultures_round2.csv')
    in_samples = metadata.loc[metadata['passage'] == 0,'sample'].values

 #   media = parent_subjects_media[-1]
   # parent_subject = parent_subjects_media[:-1]
    info, depth, freq = load_and_sort_files(species_dir, species)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]#used all for now 
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  

    depth_filtered= depth_filtering(depth,depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    s1s = []
    s2s = []
    snps_switch = []
    num_sites = []
    all_dfs = []

    for i, s1 in enumerate(freq_filtered.columns.values):
        fixed_diffs= get_fixed_diffs(freq_filtered,s1).reset_index()
        fixed_diffs['sample1'] = s1
        fixed_diffs['species_id'] = species
        all_dfs.append(fixed_diffs)

    ss_df = pd.concat(all_dfs,axis=0)
  #  ss_df['Strain Shift'] = ss_df['fixed_diffs'] > 1000
  #  if '/' in parent_subjects_media:
   #     parent_subjects_media = ''.join(parent_subjects_media.split('/'))
    ss_df.to_csv(f'{save_dir}/{species}_fixed_diffs.csv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')

#    parser.add_argument('--parent_subject_media', action = 'store', 
 #                      help = 'parent_subject-media')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species, )