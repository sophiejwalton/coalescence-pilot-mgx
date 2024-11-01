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

def get_samples_to_test():
    pass 


def get_distinguishing_snps(freq_inoculumns, thresh = .8):
    detect_df = freq_inoculumns.copy()
    detect_df = np.nan 
    detect_df[freq_inoculumns < 1- thresh] = 1 # alt allele 
    detect_df[freq_inoculumns > thresh] = 0 # ref allele 
    diffs = detect_df.diff(axis = 1)[freq_incolumns.columns.values[1]]
    diffs = diffs[diffs.abs() ==1] # both sites present 


good_inoculumns = {'AA-mBHI': 'A2-e003Coalescence-mBHI-inoculumn-redo',
 'AA-mGAM': 'A2-e003Coalescence-Inoculumn-mGAM',
 'AC/PP-mBHI': 'B3-e003Coalescence-Inoculumn-mBHI',
 'AC/PP-mGAM': 'B3-e003Coalescence-Inoculumn-mGAM',
 'AE-mBHI': 'C4-e003Coalescence-mBHI-inoculumn-redo',
 'AE-mGAM': 'C4-e003Coalescence-Inoculumn-mGAM',
 'AF-mBHI': 'D5-e003Coalescence-Inoculumn-mBHI',
 'AF-mGAM': 'D5-e003Coalescence-Inoculumn-mGAM'}

 
def get_main(species_dir, save_dir, species, parent_subjects_media):
    parent_subjects_media = parent_subjects_media.split('-')
    media = parent_subjects_media[-1]
    parent_subjects = parent_subjects_media[:-1]
    info, depth, freq = load_and_sort_files(species_dir, species)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>10.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)

    inoculumns = [good_inoculumns[f'{parent_subjects[0]}-{media}', f'{parent_subjects[1]}-{media}']]
    freq_inoculumns = freq_filtered[inoculumns]



    samples_to_test = get_samples_to_test()

    haplotype1 = get_haplotype(freq_filtered, inoculumns[0])

    haplotype2 = get_haplotype(freq_filtered, inoculumns[1])


    for i, s1 in enumerate(good_inoculumns):
        print(i, s1)
        if s1 not in good_samples.index.values:
            continue 

        for s2 in good_samples.index.values: 
            print(s2)
            if s1 == s2:
                continue 
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
    parser.add_argument('--parent_subjects', action = 'store', 
                       help = 'parent_subject-media')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species, args.parent_subjects_media)
    


       
        
    

