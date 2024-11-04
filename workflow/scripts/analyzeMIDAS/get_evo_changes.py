import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
from tqdm import tqdm
import argparse

from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
warnings.filterwarnings('ignore')

 
def get_main(subjects, save_dir, species,  main_study = True):
    changing_sites = np.zeros(len(subjects))
    good_changing_sites = np.zeros(len(subjects))
    for i, subject in enumerate(subjects):
        good_depth = pd.read_csv(f'{save_dir}/{subject}_good_depth.csv.gz', compression='gzip',).set_index('site_id')
        good_freq = pd.read_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip').set_index('site_id')
        times = np.array([int(sample.split('-')[-1]) for sample in good_freq.columns.values])
        if main_study:
            times= times[times < 75]
            samples = [f'HouseholdTransmission-Stool-{subject}-{str(time).zfill(3)}' for time in times]
            good_freq = good_freq[samples]
            good_depth = good_depth[samples]
        
        median_depth_series = good_depth.copy().replace(0, np.nan).median(skipna = True)
        
        freq_polarized_transition = get_transition_frequency_snps(good_freq.copy(), good_depth.copy())
        if len(freq_polarized_transition) == 0:
            changing_sites[i] = 0
            good_changing_sites[i] = 0
           
        depth_transition = good_depth.loc[freq_polarized_transition.index.values, :]
        freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_transition, median_depth_series)
        depth_transition_filter = depth_transition.loc[freq_transition_filter.index.values, :]
        
        changing_sites[i] = len(depth_transition)
        good_changing_sites[i] = len(depth_transition_filter)
        if main_study:
            depth_transition_filter.to_csv(f'{save_dir}/{subject}_depth_transition_main_study.csv')
            freq_transition_filter.to_csv(f'{save_dir}/{subject}_freq_transition_main_study.csv')
        else:
            depth_transition_filter.to_csv(f'{save_dir}/{subject}_depth_transition_full_study.csv')
            freq_transition_filter.to_csv(f'{save_dir}/{subject}_freq_transition_full_study.csv')

    transition_df = pd.DataFrame(data = {'Subject': subjects, 'Good Changing Sites': good_changing_sites, 
                                    'Changing sites': changing_sites})
    transition_df['Species'] = species
    if main_study:
        transition_df.to_csv(f'{save_dir}/transition_snps_df_main_study.csv')
    else:
        transition_df.to_csv(f'{save_dir}/transition_snps_df_full_study.csv')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get changing SNPS within populations')

    # add arguments

    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')
    parser.add_argument('species', action = 'store', 
                       help = 'species to perform analysis on')
    parser.add_argument('--main_study',
                    action='store_true')

    args = parser.parse_args()

    save_dir = f'{args.save_dir}/{args.species}'

    subject_fnames = glob(f'{save_dir}/*good_freq.csv.gz')
    subjects = [ fname.split('/')[-1].split('_')[0] for fname in subject_fnames]

    get_main(subjects, save_dir, args.species,  main_study = args.main_study)
    


        
    

