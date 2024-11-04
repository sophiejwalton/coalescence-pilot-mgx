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

 
def get_main(subjects, species_dir, save_dir, species):
    info, _, _ = load_and_sort_files(species_dir)
    info = info.set_index('site_id')
    changing_sites = np.zeros(len(subjects))
    good_changing_sites = np.zeros(len(subjects))

    all_strain_shift_dfs = []
    for i, subject in enumerate(subjects):
        print(subject)
        try:
            good_depth = pd.read_csv(f'{save_dir}/{subject}_good_depth.csv.gz', compression='gzip')
            good_freq = pd.read_csv(f'{save_dir}/{subject}_good_freq.csv.gz', compression='gzip')
        except:
            print('SHIT')
            continue
        if len(good_freq) == 0:
            continue 
    
        good_depth = good_depth.set_index('site_id')
        good_freq = good_freq.set_index('site_id')
        good_samples = good_freq.columns.values
        if len(good_samples) < 2:
            continue

        median_depth_series = good_depth.copy().replace(0, np.nan).median(skipna = True)
        t1s = []
        t2s = []
        snps_switch = []
        median_freq_before_switch = []
        percent_snps_at_zero = []
        s75th_before_switch = []
        s90th_before_switch = []
        focal_times = np.array([int(sample.split('-')[-1]) for sample in good_samples])
        for t1,t2 in it.combinations(focal_times, 2): 
        
            if t1 < t2:
                t_early = t1 
                t_late = t2
            else: 
                t_early = t2
                t_late = t1
            s_early = f'HouseholdTransmission-Stool-{subject}-{str(t_early).zfill(3)}'
            s_late = f'HouseholdTransmission-Stool-{subject}-{str(t_late).zfill(3)}'
            t1s.append(t_early)
            t2s.append(t_late)
            freq_small = good_freq[[s_early, s_late]].copy()
            depth_small =  good_depth[[s_early, s_late]].copy()
            freq_small = polarize_species(freq_small.copy(), s_early)
            freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
            depth_small = depth_small.loc[freq_polarized_transition.index.values, ]
            freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_small, median_depth_series[[s_early, s_late]])
            snps_switch.append(len(freq_transition_filter))

            if len(freq_transition_filter) == 0: 
                median_freq_before_switch.append(np.nan)
                percent_snps_at_zero.append(np.nan)
                s75th_before_switch.append(np.nan)
                s90th_before_switch.append(np.nan)
            else:
                freq_before_values = freq_transition_filter[s_early].values
                freq_before = np.median(freq_before_values)
                median_freq_before_switch.append(freq_before)
                print(float(freq_before))
                s75th_before_switch.append(np.percentile(freq_before_values, 75))
                s90th_before_switch.append(np.percentile(freq_before_values, 90))
                percent_zero = np.sum(freq_transition_filter[[s_early]] == 0.).values[0]/len(freq_transition_filter)
                percent_snps_at_zero.append(percent_zero)


        ss_df = pd.DataFrame(data = {'t1': t1s, 't2': t2s, 'snps_switching': snps_switch, 'median_freq': median_freq_before_switch, 
                                     '75th_Percentile': s75th_before_switch, 
                                     '90th_Percentile': s90th_before_switch, 
                                     'percent_snps_at_zero': percent_snps_at_zero })
        ss_df['Subject'] = subject
        ss_df['Species'] = species
        print(ss_df)
        all_strain_shift_dfs.append(ss_df)

    if len(all_strain_shift_dfs) == 0:
        all_strain_shift_dfs = pd.DataFrame(data = {'t1':[], 't2': [], 'snps_switching': [], 'median_freq': [], 
                                     '75th_Percentile': [], 
                                     '90th_Percentile': [], 
                                     'percent_snps_at_zero': [] })
    else:
        all_strain_shift_dfs = pd.concat(all_strain_shift_dfs)
    all_strain_shift_dfs.to_csv(f'{save_dir}/checking_pre_existing_before_shift.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('out_dir', action='store',
                    help='Outdir prefix where files are stored')
    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')
    parser.add_argument('species', action = 'store', 
                       help = 'species to perform analysis on')

#    parser.add_argument('--main_study',
 #                   action='store_true')

    args = parser.parse_args()
    species_dir = f'{args.out_dir}/{args.species}'
    save_dir = f'{args.save_dir}/{args.species}'

    subject_fnames = glob(f'{save_dir}/*good_freq.csv.gz')
    subjects = [ fname.split('/')[-1].split('_')[0] for fname in subject_fnames]
#    key_timepoint_df = pd.read_csv(args.keytimepoint_df)
    get_main(subjects, species_dir, save_dir, args.species)
    


        
    

