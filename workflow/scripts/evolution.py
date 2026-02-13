import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
#import itertools as it
#from scipy.stats import linregress
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
from track_snps_funcs import *
import warnings
warnings.filterwarnings('ignore')

def get_child_samples(inoculumn, metadata):
    subject,_,media = inoculumn.split('-')
    mesocosms = metadata.loc[metadata['parent_media'] == media,:]
    mesocosms['good']=mesocosms['parent_subjects'].transform(lambda x: subject in x)
    mesocosms = mesocosms.loc[mesocosms['good'],:]
    return mesocosms

# right now just look for mutations that are 
def get_mutations(freqs, depth, metadata, parent_dir):
    all_dfs = []
    metadata = metadata.loc[(metadata['passage']==0)+(metadata['passage']==7),:]
    metadata_ino = metadata.loc[metadata['is_inoculumn'],:]
    metadata_ino.loc[metadata_ino['parent_subjects'].isin(['AA-AA','AE-AE',
                                                                 'AF-AF']),:]
    metadata_non_zero=metadata.loc[metadata['passage']==7,:]
 #   print(metadata_non_zero',woo)
    all_mutations = []
    good_inos = np.intersect1d(freqs.columns.values, metadata_ino.index.values)
    for ino_sample in good_inos:
        inoculumn = metadata_ino.loc[ino_sample, 'inoculumn']
        metadata_good = get_child_samples(inoculumn, metadata)
        metadata_p7 = metadata_good.loc[metadata_good['passage'] >=5,:]
        good_samples = np.intersect1d(freqs.columns.values, metadata_p7.index.values)
        cand_snps = []
        for samp in good_samples:
            moving_snps = get_transition_frequency_snps(freqs[[ino_sample,samp]], depth[[ino_sample,samp]]).index.values
            if len(moving_snps) < 100:
                cand_snps = cand_snps + list(moving_snps)
        cand_snps = np.unique(cand_snps)
        good_freqs = freqs.loc[cand_snps, [ino_sample] + list(good_samples)]
        good_depth = depth.loc[cand_snps, [ino_sample] + list(good_samples)]
        good_freqs.to_csv(f'{parent_dir}/{inoculumn}_freqs_changing.csv')
        good_depth.to_csv(f'{parent_dir}/{inoculumn}_depths_changing.csv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
#    parser.add_argument('--inoculumn', action = 'store', 
 #                      help = 'inoculumn')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'

#    parent_samples, child_samples = get_parent_children(args.inoculumn)


    if not path.isdir(save_dir):
        mkdir(save_dir) 
    info, depth, freq = load_and_sort_files(species_dir, args.species)
    #print(info.columns.values)
    #print(info.index.values)
    freq = repolarize_against_reference(freq, info)
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4.csv').set_index('sample')
    
   # metadata = pd.read_csv('workflow/analysis/e003_metadata_cultures_round2_change_AA.csv')
    

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)

    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
#    late=(6,7)
 #   early = (int(args.interval[0]),int(args.interval[1]))
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)

    
   # mesos = ['A2-AA-AA-mBHI-mBHI', 'A2-AA-AA-mBHI-mGAM', 'A2-AA-AA-mGAM-mGAM',
    #   'C4-AE-AE-mBHI-mBHI', 'E2-AA-AA-mBHI-mGAM', 'E8-AA-AA-mGAM-mGAM',
     #  'G4-AE-AE-mBHI-mBHI']
   # type_mesos = ['AA-AA-mGAM-mGAM','AA-AE-mGAM-mGAM','AE-AE-mGAM-mGAM']
   # metadata = metadata.loc[metadata['type_mesocosm'].isin(type_mesos),:]
   # e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
    metadata['AC']=metadata['parent_subjects'].transform(lambda x: 'AC' in x)
    metadata= metadata.loc[~metadata['AC'],:]
    samples = np.intersect1d(metadata.index.values, freq_filtered.columns.values)
    freq_filtered = freq_filtered.loc[:,samples]
    depth_filtered = depth_filtered.loc[:,samples]
    get_mutations(freq_filtered, depth_filtered, metadata, save_dir)

       # parent2_info.to_csv(f'{save_dir}/{inoculumn}_{str(sample_init)}_{str(final_sample)}_parent2_info_shift.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

