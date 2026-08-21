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
   # metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4_good.csv')
   # in_samples = metadata.loc[metadata['passage'] == 0,'sample'].values

 #   media = parent_subjects_media[-1]
   # parent_subject = parent_subjects_media[:-1]
    info, depth, freq = load_and_sort_files(species_dir, species)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]#used all for now 
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  

    depth_filtered= depth_filtering(depth,depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    depth_filtered, freq_filtered = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
        # poll by ino 
    good_samples = np.intersect1d(freq_filtered.index.values, ['D5-e003Coalescence-Inoculumn-mGAM','A_H8_D11_AF_AF_mGAM_mGAM_2_S92',
       'B_H12_D11_AF_AF_mGAM_mGAM_4_S192',
       'C_H8_D11_AF_AF_mGAM_mGAM_6_S284', 'D11-e003Coalescence-mGAM-p5',
       'D11-e003Coalescence-mGAM-p7', 'D11-e003Coalescence-mGAM-p1_S779'])
    freq_filtered = freq_filtered[good_samples]
    freq_filtered.loc[freq_filtered['D5-e003Coalescence-Inoculumn-mGAM']>.5,:] = 1- freq_filtered.loc[freq_filtered['D5-e003Coalescence-Inoculumn-mGAM']>.5,:]
    freq_filtered = freq_filtered.loc[freq_filtered['D5-e003Coalescence-Inoculumn-mGAM']<.2,:]
    good_snvs =  freq_filtered.max(axis=1)
    good_snvs = good_snvs[good_snvs>.8].index.values
    freq_filtered.loc[good_snvs,:].to_csv('workflow/out/tryit.csv')


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