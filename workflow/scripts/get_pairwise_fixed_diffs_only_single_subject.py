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


def get_haplotype(sample_freq_filtered,sample):
    good_sites_lower = sample_freq_filtered.loc[sample_freq_filtered[sample] < .2]
    good_sites_upper = sample_freq_filtered.loc[sample_freq_filtered[sample] > .8]
    all_good_sites = pd.concat([good_sites_lower, good_sites_upper])
    all_good_sites = all_good_sites.loc[~np.isnan(all_good_sites[sample]),:]
    all_good_sites['Allele'] = 1*(all_good_sites[sample] > .5)
    return all_good_sites
 
def get_main(species_dir, save_dir, species,):
    metadata = pd.read_csv('config/e003_metadata_cultures_round2.csv')
    in_samples = metadata.loc[metadata['passage'] == 0,'sample'].values
    

 #   media = parent_subjects_media[-1]
   # parent_subject = parent_subjects_media[:-1]
    info, depth, freq = load_and_sort_files(species_dir, species)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  

    depth_filtered= depth_filtering(depth,depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    s1s = []
    s2s = []
    snps_switch = []
    num_sites = []
    #good_ins = np.intersect1d(in_samples, good_samples.index.values)
    glyc_samples = []
    for sample in good_samples.index.values:
        if 'lycero' in sample:
            glyc_samples.append(sample)
    
    good_ins = np.intersect1d(good_samples.index.values, glyc_samples+list(in_samples))
    depth_filtered = depth_filtered[good_ins]
    freq_filtered = freq_filtered[good_ins]
    for s1,s2 in it.combinations(good_ins, 2):
        print(s1,s2)
        freq_small = freq_filtered[[s1,s2]].copy()
        depth_small =  depth_filtered[[s1,s2]].copy()
        freq_small = polarize_species(freq_small.copy(), s1)   
        isna_sites= depth_small.isna().sum(axis=1)
        good_sites = isna_sites[isna_sites==0].index.values
            
        freq_small = freq_small.loc[good_sites,:]
        depth_small = depth_small.loc[good_sites,:]
        freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
        snps_switch.append(len(freq_polarized_transition))
        num_sites.append(len(good_sites))


        s1s.append(s1)
        s2s.append(s2)
    ss_df = pd.DataFrame(data = {'s1': s1s, 's2': s2s, 'fixed_diffs': snps_switch, 'num_sites': num_sites})
    ss_df['Species'] = species
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
    


       
        
    

