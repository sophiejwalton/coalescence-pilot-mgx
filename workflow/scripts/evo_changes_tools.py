import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
#from tqdm import tqdm
import argparse

from snp_analysis_tools_sherlock import *
import warnings
warnings.filterwarnings('ignore')

    
def cleanup_and_polarize(freq, median_depth_series, diversity_series, subject):
    good_samples = list(freq.columns.values)
    good_timepoints = []
    ok_timepoints = []
    times = []

    for sample in good_samples:
        if median_depth_series[sample]>= 20:
            if diversity_series[sample] < 1e-3:
                good_timepoints.append(int(sample.split('-')[-1]))
            else:
                ok_timepoints.append(int(sample.split('-')[-1]))
        else: 
          #  print('cool')
            times.append(int(sample.split('-')[-1]))
  #  print(len(good_samples), len(good_timepoints), len(ok_timepoints), len(times))
    if len(good_timepoints) > 0:
        first_good_time = np.min(np.array(good_timepoints))
    elif len(ok_timepoints)> 0:
        first_good_time = np.min(np.array(ok_timepoints))
    else: 
        first_good_time = np.min(times)
        
   # print('YAY', good_samples)
    sample_to_polarize = f'HouseholdTransmission-Stool-{subject}-{str(first_good_time).zfill(3)}'
   # print(sample_to_polarize)
    freq_polarized = polarize_species(freq[good_samples].copy(), sample_to_polarize)
    return freq_polarized, good_samples


def filter_sites_across_samples(good_depth, good_freq):
    
    good_samples = good_depth.columns
    counts = good_depth.count(axis = 1)
    
    mintimes = round(len(good_samples)*.8)
    passing_sites = counts[counts > mintimes]
       # print(len(good_freq))
    good_freq = good_freq.loc[passing_sites.index.values, :]
    good_depth = good_depth.loc[passing_sites.index.values, :]
       # print(len(good_freq))
    return good_depth, good_freq

def get_transition_frequency_snps(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, .5) # so do not help o

    freq_pass_2 = freq_polarized[(temp_freq<0.2).any(axis=1)]
    depth_pass_2 = depth_filtered[(temp_freq<0.2).any(axis=1)]
  #  print(len(freq_pass_2))
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
        
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, .5)

    freq_polarized_transition= freq_pass_2[(temp_freq > 0.8).any(axis=1)]
    #print(freq_polarized_transition)
    return freq_polarized_transition


def filter_transition_frequency(moving_snps,  depth_transition, median_depth_series):
    good_sites = []
    for site in moving_snps.index.values:
        good_changing_site = False 
      
        
        coverage_ratio = depth_transition.loc[site, :]/median_depth_series[depth_transition.columns]
        coverage_ratio = coverage_ratio
        site_freq = moving_snps.loc[site, :]
        low_samples  = site_freq[site_freq < .2]
        high_samples = site_freq[site_freq > .8]
        high_samples_coverage = coverage_ratio[site_freq > .8]
        low_samples_coverage = coverage_ratio[site_freq < .2]
        if high_samples_coverage.max()/low_samples_coverage.min() < 2:
            good_changing_site = True 

        else:

            for sample in low_samples.index:
                good_highs = high_samples_coverage[high_samples_coverage < coverage_ratio[sample]*2]

                if len(good_highs) > 0:
                    good_changing_site = True

                    break
        if good_changing_site:
                good_sites.append(site)
    
    freq_filtered_transition = moving_snps.loc[good_sites, :]
    
    return  freq_filtered_transition 


  