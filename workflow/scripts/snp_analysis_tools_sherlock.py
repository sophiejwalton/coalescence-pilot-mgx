import pandas as pd
import numpy as np
from os import path
#import lz4.frame

from collections import Counter


from glob import glob
#from tqdm import tqdm

def load_and_sort_files(species_dir,species):
    """
    Takes in snp directory, returns dataframes of various MIDAS2 outputs: info, depth, freq
    """
    info = pd.read_csv(species_dir + f"/{species}.snps_info.tsv.gz", sep = "\t", compression='gzip', low_memory=False).set_index('site_id')
    depth = pd.read_csv(species_dir + f"/{species}.snps_depth.tsv.gz", compression='gzip',sep = "\t", ).set_index('site_id')
    freq = pd.read_csv(species_dir + f"/{species}.snps_freqs.tsv.gz",compression='gzip', sep = "\t", ).set_index('site_id')

    return info, depth, freq



def depth_filtering(depth):
    depth_no_site = depth.copy()
    med = depth_no_site.replace(0, np.nan).median(axis = 0,skipna=True)
    depth_masked_1 = depth_no_site.mask((depth_no_site > 2.5 * med),axis = 0)
    depth_masked = depth_masked_1.mask((depth_masked_1 < med / 2.5),axis = 0)
    depth_masked_absolute = depth_masked.mask(depth_masked < 5)
    return depth_masked_absolute

def freq_masked(freq, depth_filtered):
    depth_filtered_na = depth_filtered*depth_filtered.isna() + 1. 
    return freq.copy()*depth_filtered_na 


def get_diversity_series(freq_filtered, thresh=.2):
    temp_freq = freq_filtered[freq_filtered.notna()].replace(np.nan, .5)
    less_than_thresh = temp_freq < thresh
    less_than_thresh = less_than_thresh.sum()
    greater_than_thresh = temp_freq > 1-thresh
    greater_than_thresh = greater_than_thresh.sum()
    all_nonint_sites = greater_than_thresh + less_than_thresh
    # number of intermediate frequency sites = number total sites with non nan depth - non intermediate frequency sites 
    num_int_sites = freq_filtered.count() - all_nonint_sites
    diversity = num_int_sites/freq_filtered.count() 
    return num_int_sites, diversity
    
def polarize_species(freq, sample):
    """
    Polarizes all time points to a single sample. If the frequency in the time point to polarize is >0.5,
    all other timepoints at that site are converted to 1-freq.
    Inputted frequency should not have site_id as a column, but as the index.
    """

    samples = list(freq.columns)

    if 'site_id' not in freq.columns:
        freq.index.name = 'site_id'
        # Get site id as a column
        freq.reset_index(inplace = True)

    # Turn into a melted df
    freq_melted = pd.melt(freq, id_vars=['site_id'], value_name = 'freq', var_name = 'sample')
#    print(freq_melted.head())
    freq_melted_polarized = freq_melted.copy()
    sites_to_flip = freq.loc[freq[sample] > 0.5]['site_id']

    freq_melted_polarized['freq'].where(~freq_melted['site_id'].isin(sites_to_flip), 1 - freq_melted['freq'], inplace = True)

    # Return to rectangle dataframe
    freq_polarized = freq_melted_polarized.set_index(['site_id', 'sample'])['freq'].unstack()
  #  print(freq_polarized.head())
    if 'site_id'  in samples:
        samples.remove('site_id')

    # Reorder columns
  #  print(freq_polarized.columns)
#    freq_polarized = freq_polarized[samples]
    
    return freq_polarized


def get_qp_sites(freq, depth_filtered, genome_length, sites_considered):
    temp_freq = freq[depth_filtered.notna()].replace(np.nan, .5)

    less_than_20 = temp_freq < .2
    less_than_20 = less_than_20.sum()

    greater_than_80 = temp_freq > .8
    greater_than_80 = greater_than_80.sum()

    all_nonint_sites = greater_than_80 + less_than_20 
    # number of intermediate frequency sites = number total sites with non nan depth - non intermediate frequency sites 
    num_int_sites = depth_filtered.count() - all_nonint_sites

    # number of bad sites = number sites considered - number good sites 
    bad_sites = sites_considered - depth_filtered.count()

    # total sites = genome length - bad sites 
    total_sites = -bad_sites + genome_length
    qp_sites = num_int_sites/total_sites
    return qp_sites, num_int_sites 
  



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


def repolarize_against_reference(freq, info):
    '''
    Repolarize so that all allele frequencies are the frequency of the alternative allele 
    '''
    info2 = info.reset_index().copy()
    info2['site_id_adjust'] = info2['site_id'] 
    info2 = info2.set_index('site_id_adjust')
    repolarize = info2['ref_allele'] == info2['minor_allele']

    repolarized = repolarize.copy()
    repolarized[:] = False
    
    repolarized_index = repolarize.loc[repolarize == True].index 
   
    repolarized[repolarized_index] = True
    
    freq_polarized = freq.copy().set_index('site_id')
    
    freq_polarized.loc[repolarized,: ] = 1 - freq_polarized.loc[repolarized,: ]
    return freq_polarized,repolarized, repolarize

    

def get_intermediate_frequency_snps(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, -1)

    freq_pass_1 = freq_polarized[~(temp_freq<0.0).any(axis=1)]
    depth_pass_1 = depth_filtered[~(temp_freq<0.0).any(axis=1)]
   
    freq_pass_2 = freq_pass_1[~(freq_pass_1<0.2).all(axis=1)]
    depth_pass_2 = depth_pass_1[~(freq_pass_1<0.2).all(axis=1)]
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
  
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, 1.1)
    freq_polarized_plotting = freq_pass_2[~(temp_freq > 0.8).all(axis=1)]
    return freq_polarized_plotting


def get_intermediate_frequency_snps_v2(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, -1)

    freq_pass_2 = freq_polarized[~(temp_freq<0.2).all(axis=1)]
    depth_pass_2 = depth_filtered[~(temp_freq<0.2).all(axis=1)]
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
  
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, 1.1)

    freq_polarized_plotting = freq_pass_2[~(temp_freq > 0.8).all(axis=1)]
    return freq_polarized_plotting



    
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


  
