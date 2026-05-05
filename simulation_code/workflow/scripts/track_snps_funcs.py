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


def get_distinguishing_snps(freq_inoculumns, thresh = .8):
    # find snps that are greater than .8 (aka alternative alleles > .8)
    detect_df = freq_inoculumns >= thresh
    detect_df_opp = freq_inoculumns <= 1 - thresh
    detect_df['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df_opp['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df = detect_df.loc[detect_df['site_present']==0,:]
    detect_df_opp = detect_df_opp.loc[detect_df_opp['site_present']==0,:]
    detect_df['diff1'] = detect_df[freq_inoculumns.columns.values[0]].astype(int) + detect_df_opp[freq_inoculumns.columns.values[1]].astype(int)
    detect_df['diff2'] = detect_df[freq_inoculumns.columns.values[1]].astype(int) + detect_df_opp[freq_inoculumns.columns.values[0]].astype(int)
    return detect_df['diff1'], detect_df['diff2']

def filter_distinguishing_snps(freq_children, parent_snps, thresh = .5, sample_thresh=1.):
    med = freq_children.loc[parent_snps,:].median(axis = 0) 
    diff_from_med = freq_children - med
   # print(med)
    #print('PARTY')
    #print(diff_from_med)
    freq_masked = freq_children.mask((diff_from_med.abs() > thresh),axis=0)

    min_samples = round(len(freq_children.columns.values)*sample_thresh)
    passing_sites = freq_masked.count(axis = 1)
    passing_sites = passing_sites[passing_sites>= min_samples].index.values
    parent_snps_good = np.intersect1d(passing_sites, parent_snps)
    return parent_snps_good



def get_frequency_parent(freq_children, parent_snps):
    median_freq = freq_children.loc[parent_snps,:].median(axis = 0)
    return median_freq 

def get_quantile_parent(freq_children, parent_snps, q= 97.5):
    median_freq = freq_children.loc[parent_snps,:].quantile(q, axis = 0)
    return median_freq 

  # DataFrame.quantile(q=0.5, axis=0, numeric_only=False, interpolation='linear', method='single')

def get_frequency_parent_avg(freq_children, parent_snps, freq_thresh = 1.5):
    #med = freq_children.loc[parent_snps,:].median(axis = 0)
    #freq_masked = freq_children.mask((freq_children > freq_thresh * med),axis = 0)
   # freq_masked = freq_masked .mask((freq_masked  < med / freq_thresh),axis = 0)
    freq_avg = (freq_children.loc[parent_snps,:].sum(axis = 0))/(len(parent_snps))
    return freq_avg

def get_frac_zero(freq_children, parent_snps):
    frac_zero = (freq_children.loc[parent_snps,:] == 0).sum(axis=0)/len(parent_snps)
    return frac_zero 

def fix_zeros(freq_parent,  depth_parent, freq_children, parent_snps):
    freq_parent_fix = freq_parent.copy()
    frac_zero_parent= get_frac_zero(freq_children, parent_snps)
   # print(frac_zero_parent,'w0')
    freq_parent_fix[freq_parent==0] = -np.log(frac_zero_parent[freq_parent==0])/depth_parent[freq_parent==0]
    #print(np.log(1-frac_zero_parent[freq_parent==0]))
    freq_parent_pol = 1-freq_parent.copy()
    frac_one_parent= get_frac_zero(1-freq_children.copy(), parent_snps)
    freq_parent_fix[freq_parent_pol==0] = 1 + np.log(frac_one_parent[freq_parent_pol==0])/depth_parent[freq_parent_pol ==0]
    return freq_parent_fix

def fix_zeros_bs(sample_meds, depth, bs_samples,n_snps):
    freq_parent_fix = sample_meds.copy()
    frac_zero_parent= np.sum(bs_samples == 0, axis=0)/n_snps
    print(len(frac_zero_parent), np.max(frac_zero_parent),n_snps, len(freq_parent_fix))
    freq_parent_fix[freq_parent_fix == 0] = -np.log(frac_zero_parent[freq_parent_fix==0])/depth 
    frac_one_parent= np.sum(bs_samples == 1, axis=0)/n_snps
    print(np.max(frac_one_parent),'pot')
    freq_parent_fix[freq_parent_fix==1] = 1 + np.log(frac_one_parent[freq_parent_fix==1])/depth
    return freq_parent_fix

def get_freq_est(snps,depth):

    med = np.nanmedian(snps)
    n_snps = len(snps)- np.sum(np.isnan(snps))
    if med == 0:
        frac_zero = np.sum(snps == 0)/n_snps
        return -np.log(frac_zero)/depth
    elif med == 1:
        frac_one = np.sum(snps == 1)/n_snps
        return 1+ np.log(frac_one)/depth
    return med 


def get_bootstrap_parent(freq_children, depth_med, parent_snps, n_bootstraps = 1000):
    #med = freq_children.loc[parent_snps,:].median(axis = 0)
    #freq_masked = freq_children.mask((freq_children > freq_thresh * med),axis = 0)
   # freq_masked = freq_masked .mask((freq_masked  < med / freq_thresh),axis = 0)
    snps = freq_children.loc[parent_snps,:]
    boot_med = []
    boot_low = []
    boot_high = []
    act_med = []
    samples_good = []
    for sample in snps.columns.values:
        snps_sample= snps[sample].values
        
        snps_sample = snps_sample[~np.isnan(snps_sample)]
        med_og = np.median(snps_sample)
        #print(snps_sample)
        n_snps = len(parent_snps)
        print(med_og)
     #   print(np.max(snps_sample))
        if med_og == 0:
            med_og = -np.log(np.sum(snps_sample == 0)/n_snps)/depth_med[sample]
            print('zero', med_og) 
        elif med_og == 1:
            med_og = 1+np.log(np.sum(snps_sample == 1)/n_snps)/depth_med[sample] 
            print('one', med_og)
    #    act_med.append(med_og)
        if len(snps_sample) == 0: continue 
        bs_samples = np.random.choice(snps_sample, size = (len(snps), n_bootstraps))
        samples_meds = np.median(bs_samples,axis = 0)
#species_list = [102279]
#inoculumns = ['AA-AE-mGAM', 'AA-AF-mGAM'
        print(np.sum(samples_meds==1),np.sum(samples_meds==0), 'sums')  
        samples_meds_fix = fix_zeros_bs(samples_meds, depth_med[sample], bs_samples,n_snps)
        samples_good.append(sample)
        boot_med.append(np.median(samples_meds_fix))
        boot_low.append(np.percentile(samples_meds_fix, q =2.5))
        act_med.append(med_og)
        boot_high.append(np.percentile(samples_meds_fix, q = 97.5))
    return pd.DataFrame(data = {'sample': samples_good, 'boot_med': boot_med, 'boot_low': boot_low, 'boot_high': boot_high,
                                    'actual_med': act_med})


def get_parent_children(inoculumn, metadata):
#    metadata = pd.read_csv('config/e003_coal_metadata_full.csv')
    
    child_samples = list(metadata.loc[metadata['inoculumn'] == inoculumn, 'sample'].values)
    
    parent_subjects = inoculumn.split('-')[:-1]
    parent_media = inoculumn.split('-')[-1]
    ins = metadata.loc[metadata['is_inoculumn'],:]
    ins = ins.loc[ins['parent_media'] == parent_media,:]
    print(parent_subjects)
    ins1 = ins.loc[ins['parent_subjects'] == parent_subjects[0] + '-' +  parent_subjects[0],:]

    ins2 = ins.loc[ins['parent_subjects'] == parent_subjects[1] + '-' +  parent_subjects[1],:]
    if (len(ins1) ==0) or (len(ins2) ==0): 
        return np.nan, np.nan
    parent_samples = [ins1['sample'].values[0], ins2['sample'].values[0]]
    media = inoculumn.split('-')[-1]
    in_parent1 = f'{parent_subjects[0]}-{parent_subjects[0]}-{media}'
    in_parent2 = f'{parent_subjects[1]}-{parent_subjects[1]}-{media}'
  #  print(in_parent1, in_parent2)
    child_samples_parent1_ss =  list(metadata.loc[metadata['inoculumn'] == in_parent1, 'sample'].values)
    child_samples_parent2_ss =  list(metadata.loc[metadata['inoculumn'] == in_parent2, 'sample'].values)
    #print(child_samples + child_samples_parent1_ss + child_samples_parent2_ss )
    child_samples_all = child_samples # + child_samples_parent1_ss + child_samples_parent2_ss 
    return parent_samples, child_samples_all
     

def get_depth_parent(depth_children, parent_snps):
    sum_freq = depth_children.loc[parent_snps,:].sum(axis = 0,skipna=True)

    return sum_freq 


def get_count(freq_children, parent_snps):
    count_parent = freq_children.loc[parent_snps,:].count(axis=0)
    return count_parent
