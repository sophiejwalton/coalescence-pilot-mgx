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

 
def get_main(subject_species_df, subject, save_dir, key_timepoint_df):
    subject_species_df['t1_t2'] = subject_species_df['t1'].astype(str) + '_' + subject_species_df['t2'].astype(str)
    good_depth = pd.read_csv(f'{save_dir}/{subject}_good_depth.csv.gz', compression='gzip').set_index('site_id')
    median_depth_series = good_depth.copy().replace(0, np.nan).median(skipna = True)
    good_freq = pd.read_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip').set_index('site_id')
    good_samples_in= (key_timepoint_df.loc[key_timepoint_df['Subject'] == subject, 'Sample'].values)
    good_samples = [f'HouseholdTransmission-Stool-{sample}' for sample in good_samples_in]
  
    focal_samples = []
    for sample in good_samples:
        if sample in good_freq.columns.values: 
            focal_samples.append(sample)
    focal_times = np.array([int(sample.split('-')[-1]) for sample in focal_samples])
    median_max_check = []
    samples_before_num = []
    good_freq = good_freq[focal_samples]
    good_depth = good_depth[focal_samples]
    for i, t1_t2_pair in enumerate(subject_species_df['t1_t2']):
        small_df = subject_species_df.loc[subject_species_df['t1_t2'] == t1_t2_pair, :]
        t1 = small_df['t1'].values[0]
        t2 = small_df['t2'].values[0]
        s_early = f'HouseholdTransmission-Stool-{subject}-{str(t1).zfill(3)}'
        s_late = f'HouseholdTransmission-Stool-{subject}-{str(t2).zfill(3)}'
        freq_small = good_freq[[s_early, s_late]].copy()
        
        freq_polarized = polarize_species(good_freq.copy(), s_early)
        freq_small = freq_polarized[[s_early, s_late]].copy()
        depth_small =  good_depth[[s_early, s_late]].copy()

        freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
        freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_small, median_depth_series[[s_early, s_late]])

       # print('first',freq_polarized.median())
        print('transition', freq_transition_filter.median())
        print('transition', freq_transition_filter.count())
        times_before = focal_times[focal_times < t1]
  
        samples_before = [f'HouseholdTransmission-Stool-{subject}-{str(time).zfill(3)}' for time in times_before]
        if len(times_before) == 0:
            median_max_check.append(np.nan)
            samples_before_num.append(0)
        else:
            freq_transition_filter_all = freq_polarized.loc[freq_transition_filter.index.values, samples_before]
            freq_before_values_medians = freq_transition_filter_all.median(skipna=True).values
            print('before', freq_transition_filter_all.median(skipna=True))
            print('before', freq_transition_filter_all.count())
           # print('before', freq_transition_filter_all.index.values)
            freq_before_values_median_max = np.max(freq_before_values_medians)
            median_max_check.append(freq_before_values_median_max)
            samples_before_num.append(len(times_before))
    subject_species_df['median_max_check'] = median_max_check
    subject_species_df['number of samples checked'] = samples_before_num
    def get_check_value(x):
        if np.isnan(x):
            return np.nan
        if x == 0:
            return False
        if x > 0:
            return True 
    subject_species_df['evidence for preexisting'] =  subject_species_df['median_max_check'].transform(get_check_value)
    return subject_species_df

    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')
    parser.add_argument('population_df', action = 'store', 
                       help = 'file with strain shift and preexiting info')
    parser.add_argument('keytimepoint_df', action = 'store', 
                       help = 'keytimepoint_df_loc')

#    parser.add_argument('--main_study',
 #                   action='store_true')

    args = parser.parse_args()
    population_df = pd.read_csv(args.population_df)
    population_df_nonpre = population_df.loc[population_df['median_freq'] == 0.0, :]
    species_list = population_df_nonpre['Species'].unique()
    key_timepoint_df = pd.read_csv(args.keytimepoint_df)
    all_annotated_dfs = []

    for species in species_list:
        species_df = population_df_nonpre.loc[population_df_nonpre['Species'] == species, :]
        subjects = species_df['Subject'].unique()
        for subject in subjects:
            subject_species_df = species_df.loc[species_df['Subject'] == subject, :]
            subject_species_df_annotated = get_main(subject_species_df, subject, f'{args.save_dir}/{species}', key_timepoint_df)
            all_annotated_dfs.append(subject_species_df_annotated)
    all_anotated_dfs = pd.concat(all_annotated_dfs)
    all_anotated_dfs.to_csv(f'workflow/reports/preexisting_recheck_focal.csv')
    


        
    

