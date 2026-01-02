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


# right now just look for mutations that are 
def get_mutations(freqs, depth, metadata):
    all_dfs = []
    metadata = metadata.loc[(metadata['passage']==0)+(metadata['passage']==7),:]
    metadata_non_zero=metadata.loc[metadata['passage']==7,:]
    all_mutations = []
    for mesocosm in metadata['mesocosm'].unique():
        meso_df = metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, :]
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        freqs_good = freqs[[inoculumn]+ list(meso_df['sample'].values)]
        freqs_polarize = freqs_good.copy()
        freqs_polarize.loc[freqs_polarize[inoculumn]>.8,:] = 1- freqs_polarize.loc[freqs_polarize[inoculumn]>.8,:]
        print(freqs_polarize)
        samples7 =  meso_df.loc[meso_df['passage']>=7,'sample']
        mutations = freqs_polarize[samples7]>.8
        mutations = mutations[mutations].index.values
        all_mutations = all_mutations + list(mutations)

    all_mutations = np.unique(all_mutations) 
    return freqs.loc[all_mutations,:], depth.loc[all_mutations,:]



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
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4.csv')
    
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

    
    mesos = ['A2-AA-AA-mBHI-mBHI', 'A2-AA-AA-mBHI-mGAM', 'A2-AA-AA-mGAM-mGAM',
       'C4-AE-AE-mBHI-mBHI', 'E2-AA-AA-mBHI-mGAM', 'E8-AA-AA-mGAM-mGAM',
       'G4-AE-AE-mBHI-mBHI']
    metadata = metadata.loc[metadata['mesocosm'].isin(mesos),:]
    samples = np.intersect1d(metadata['sample'].values, freq_filtered.columns.values)
    freq_filtered = freq_filtered.loc[:,samples]
    depth_filtered = depth_filtered.loc[:,samples]
    freq_small, depth_small = get_mutations(freq_filtered, depth_filtered, metadata)

    freq_small.to_csv(f'{save_dir}/freq_small_evo.csv')
    depth_small.to_csv(f'{save_dir}/depth_small_evo.csv')
       # parent2_info.to_csv(f'{save_dir}/{inoculumn}_{str(sample_init)}_{str(final_sample)}_parent2_info_shift.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

