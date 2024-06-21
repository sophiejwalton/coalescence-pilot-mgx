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

 
def get_main(species_dir, save_dir, species):
    info, depth, freq = load_and_sort_files(species_dir, species)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>10.]
    depth = depth[good_samples.index.values]
    print('why')
    freq = freq[good_samples.index.values]  
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)

    s1s = []
    s2s = []
    snps_switch = []
    print('yay', good_samples.index.values)
    inoculumns = []
    good_inoculumns = ['A2-e003Coalescence-Inoculumn-mBHI',
       'A2-e003Coalescence-Inoculumn-mGAM',
       'A2-e003Coalescence-mBHI-inoculumn-redo',
       'B3-e003Coalescence-Inoculumn-mBHI',
       'B3-e003Coalescence-Inoculumn-mGAM',
       'C4-e003Coalescence-Inoculumn-mGAM',
       'C4-e003Coalescence-mBHI-inoculumn-redo',
       'D5-e003Coalescence-Inoculumn-mBHI',
       'D5-e003Coalescence-Inoculumn-mGAM']
    for s in good_samples.index.values:
        if s in good_inoculumns:
            inoculumns.append(s)

    for i, s1 in enumerate(inoculumns):
        print(i, s1)

        for s2 in good_samples.index.values:
            freq_small = freq_filtered[[s1,s2]].copy()
            depth_small =  depth_filtered[[s1,s2]].copy()
            freq_small = polarize_species(freq_small.copy(), s1)
            freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
            depth_small = depth_small.loc[freq_polarized_transition.index.values, ]
            freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_small, med_nonzero_depth[[s1,s2]])
            snps_switch.append(len(freq_transition_filter))
            s1s.append(s1)
            s2s.append(s2)
    ss_df = pd.DataFrame(data = {'s1': s1s, 's2': s2s, 'fixed_diffs': snps_switch, })
    ss_df['Species'] = species
    ss_df['Strain Shift'] = ss_df['fixed_diffs'] > 1000
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
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species)
    


       
        
    

