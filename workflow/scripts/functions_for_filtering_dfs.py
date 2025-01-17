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


def depth_filtering(depth, depth_thresh=2.5):
    '''
    marks which sites fail filtering in each sample based on relative coverage 
    '''
    depth_no_site = depth.copy()
    med = depth_no_site.replace(0, np.nan).median(axis = 0,skipna=True) # get median nonzero depth- important 
    depth_masked_1 = depth_no_site.mask((depth_no_site > depth_thresh * med),axis = 0)
    depth_masked = depth_masked_1.mask((depth_masked_1 < med / depth_thresh),axis = 0)
    depth_masked_absolute = depth_masked.mask(depth_masked < 0) #DO WE ABSOLUTE MAX THIS 
    return depth_masked_absolute



def freq_masked(freq, depth_filtered):
    '''
    using the depth df as a marker of which sites pass filtering, sets the value of the frequency df to nan 
    for these sites in each sample
    '''
    
    depth_filtered_na = depth_filtered*depth_filtered.isna() + 1. 
    return freq.copy()*depth_filtered_na



def filter_sites_across_samples(depth_filtered, freq_filtered, thresh = .8):
    '''
    removes sites that are missing from too many samples
    '''
    
    good_samples = depth_filtered.columns
    counts = depth_filtered.count(axis = 1)
    
    mintimes = round(len(good_samples)*thresh) # has to be fully present
    passing_sites = counts[counts > mintimes]
       # print(len(good_freq))
    good_freq = freq_filtered.loc[passing_sites.index.values, :]
    good_depth = depth_filtered.loc[passing_sites.index.values, :]
       # print(len(good_freq))
    return good_depth, good_freq 

def repolarize_against_reference(freq, info):
    '''
    Repolarize so that all allele frequencies are the frequency of the alternative allele 
    '''
    #info2 = info.copy().reset_index().set_index('site_id')
  #  freq_polarized = freq.copy().set_index('site_id')
    info2 = info.reset_index()
    info2['ref_allele'] = info2['site_id'].transform(lambda x: x.split('|')[-1])
    info2 = info2.set_index('site_id')
    freq.loc[info2['ref_allele'] == info2['minor_allele'],: ] = 1 - freq.loc[info2['ref_allele'] == info2['minor_allele'],: ]

    return freq

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

