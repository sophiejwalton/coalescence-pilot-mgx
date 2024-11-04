import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
from tqdm import tqdm
import argparse
import itertools as it
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
warnings.filterwarnings('ignore')

 
def get_main(subjects, species_dir, save_dir, species, key_timepoint_df):
    info, _, _ = load_and_sort_files(species_dir)
    info = info.set_index('site_id')
    changing_sites = np.zeros(len(subjects))
    good_changing_sites = np.zeros(len(subjects))
    all_strain_shift_dfs = []
    for i, subject in enumerate(subjects):
        good_depth = pd.read_csv(f'{save_dir}/{subject}_good_depth.csv.gz', compression='gzip',).set_index('site_id')
        good_freq = pd.read_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip').set_index('site_id')
       # times = np.array([int(sample.split('-')[-1]) for sample in good_freq.columns.values])
        good_samples = (key_timepoint_df.loc[key_timepoint_df['Subject'] == subject, 'Sample'].values)
       # print(focal_samples)
        good_samples = [f'HouseholdTransmission-Stool-{sample}' for sample in good_samples]
        focal_samples = []
        for sample in good_samples:
            if sample in good_freq.columns.values: 
                focal_samples.append(sample)
        if len(focal_samples) < 2:
            continue
        print(focal_samples)
        print(good_freq.columns.values)
        good_freq = good_freq[focal_samples]
        good_depth = good_depth[focal_samples]
        median_depth_series = good_depth.copy().replace(0, np.nan).median(skipna = True)
        t1s = []
        t2s = []
        snps_switch = []
        focal_times = np.array([int(sample.split('-')[-1]) for sample in focal_samples])
        for t1,t2 in it.combinations(focal_times, 2): 
            t1s.append(t1)
            t2s.append(t2)
            s1 = f'HouseholdTransmission-Stool-{subject}-{str(t1).zfill(3)}'
            s2 = f'HouseholdTransmission-Stool-{subject}-{str(t2).zfill(3)}'
            freq_small = good_freq[[s1, s2]].copy()
            depth_small =  good_depth[[s1, s2]].copy()
            freq_small = polarize_species(freq_small.copy(), s1)
            freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
            depth_small = depth_small.loc[freq_polarized_transition.index.values, ]
            freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_small, median_depth_series[[s1, s2]])
            snps_switch.append(len(freq_transition_filter))


        ss_df = pd.DataFrame(data = {'t1': t1s, 't2': t2s, 'snps_switching': snps_switch})
        ss_df['Subject'] = subject
        ss_df['Species'] = species
        print(ss_df)
        all_strain_shift_dfs.append(ss_df)

    if len(all_strain_shift_dfs) == 0:
        all_strain_shift_dfs = pd.DataFrame(data = {'t1': [], 't2': [], 'snps_switching': []})
    else:
        all_strain_shift_dfs = pd.concat(all_strain_shift_dfs)
    all_strain_shift_dfs.to_csv(f'{save_dir}/strain_shifts_df_focal_study.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('out_dir', action='store',
                    help='Outdir prefix where files are stored')
    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')
    parser.add_argument('species', action = 'store', 
                       help = 'species to perform analysis on')
    parser.add_argument('keytimepoint_df', action = 'store', 
                       help = 'keytimepoint_df_loc')
#    parser.add_argument('--main_study',
 #                   action='store_true')

    args = parser.parse_args()
    species_dir = f'{args.out_dir}/{args.species}'
    save_dir = f'{args.save_dir}/{args.species}'

    subject_fnames = glob(f'{save_dir}/*good_freq.csv.gz')
    subjects = [ fname.split('/')[-1].split('_')[0] for fname in subject_fnames]
    key_timepoint_df = pd.read_csv(args.keytimepoint_df)
    get_main(subjects, species_dir, save_dir, args.species, key_timepoint_df)
    


        
    

