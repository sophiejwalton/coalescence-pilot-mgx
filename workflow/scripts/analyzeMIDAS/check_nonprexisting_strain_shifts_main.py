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

 
def get_main(subject_species_df, subject, save_dir, main_study, katherine):
    subject_species_df['t1_t2'] = subject_species_df['t1'].astype(str) + '_' + subject_species_df['t2'].astype(str)
    good_depth = pd.read_csv(f'{save_dir}/{subject}_good_depth.csv.gz', compression='gzip').set_index('site_id')
    median_depth_series = good_depth.copy().replace(0, np.nan).median(skipna = True)
    good_freq = pd.read_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip').set_index('site_id')
#    good_samples = [f'HouseholdTransmission-Stool-{sample}' for sample in good_samples_in]
    good_samples_all = good_freq.columns.values
  
    good_times_all = np.array([int(sample.split('-')[-1]) for sample in good_samples_all])
    if main_study:
        good_times = good_times_all[good_times_all <= 75] 
    else:
        good_times = good_times_all
    good_samples = [ f'HouseholdTransmission-Stool-{subject}-{str(time).zfill(3)}' for time in good_times]
    median_max_check = []
    samples_before_num = []
    good_freq = good_freq[good_samples]
    good_depth = good_depth[good_samples]
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
      #  print('transition', freq_transition_filter.index.values)
        times_before = good_times[good_times < t1]
  
        samples_before = [f'HouseholdTransmission-Stool-{subject}-{str(time).zfill(3)}' for time in times_before]
        if len(times_before) == 0:
            median_max_check.append(np.nan)
            samples_before_num.append(0)
        else:
            freq_transition_filter_all = freq_polarized.loc[freq_transition_filter.index.values, samples_before]
            freq_before_values_medians = freq_transition_filter_all.median().values
            print('before', freq_transition_filter_all.median())
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

    parser.add_argument('--main_study',
                   action='store_true')

    parser.add_argument('--katherine',
                   action='store_true')

    args = parser.parse_args()
    population_df = pd.read_csv(args.population_df)
    population_df_nonpre = population_df.loc[population_df['median_freq'] == 0.0, :]
    population_df_nonpre['population'] = population_df_nonpre['Species'] + ' in ' + population_df_nonpre['Subject']
    if args.katherine:
        population_df_katherine = pd.read_csv('workflow/reports/strainTurnover-populations.txt', delimiter = ' ')
        population_df_katherine['population'] = population_df_katherine['species_id'] + ' in ' + population_df_katherine['subject']
    population_df_nonpre = population_df_nonpre.loc[population_df_nonpre['population'].isin(population_df_katherine['population'].unique()),:]
    species_list = population_df_nonpre['Species'].unique()
    big_df = []
    for species in species_list:
        species_df = population_df_nonpre.loc[population_df_nonpre['Species'] == species, :]
        subjects = species_df['Subject'].unique()
        all_annotated_dfs = []
        for subject in subjects:
            subject_species_df = species_df.loc[species_df['Subject'] == subject, :]
            subject_species_df_annotated = get_main(subject_species_df, subject, f'{args.save_dir}/{species}', args.main_study, args.katherine)
            all_annotated_dfs.append(subject_species_df_annotated)
            big_df.append(subject_species_df_annotated)
        all_anotated_dfs = pd.concat(all_annotated_dfs)
        if args.main_study:
            all_anotated_dfs.to_csv(f'workflow/reports/{species}_preexisting_recheck_main_study.csv')
        elif args.katherine:
            all_anotated_dfs.to_csv(f'workflow/reports/{species}_preexisting_recheck_katherine.csv')
        else:
            all_anotated_dfs.to_csv(f'workflow/reports/{species}_preexisting_recheck_full_study.csv')
    big_df = pd.concat(big_df)
    big_df.to_csv('workflow/reports/preexisting_recheck_katherine_all.csv')

        


            
    

