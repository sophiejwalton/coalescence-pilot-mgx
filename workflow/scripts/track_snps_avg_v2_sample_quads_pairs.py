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
    
    detect_df['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df = detect_df.loc[detect_df['site_present']==0,:]
    print(detect_df)
#    print(detect_df[freq_inoculumns.columns.values[0]])
    detect_df['diff'] = detect_df[freq_inoculumns.columns.values[0]].astype(int) - detect_df[freq_inoculumns.columns.values[1]].astype(int)
    return detect_df['diff']

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


def get_freq_est(snps,depth)
    med = np.nanmedian(snps)
    if med == 0:
        frac_zero = np.sum(snps == 0)/len(snps)
        return -np.log(frac_zero)/depth
    elif med == 1:
        frac_one = np.sum(snps == 1)/len(snps)
        return 1+ np.log(frac_one)/depth
    return med 


def get_bootstrap_sel_coeffs(metadata, freq_children, depth_med, parent_snps, n_bootstraps = 1000):
    #med = freq_children.loc[parent_snps,:].median(axis = 0)
    #freq_masked = freq_children.mask((freq_children > freq_thresh * med),axis = 0)
   # freq_masked = freq_masked .mask((freq_masked  < med / freq_thresh),axis = 0)
    snps = freq_children.loc[parent_snps,:]
    print(snps.columns.values)
    boot_med = []
    boot_low = []
    boot_high = []
    act_med = []
    samples_good = []
    metadata = metadata.loc[metadata['sample'].isin(snps.columns.values),:]
    metadata_non_zero = metadata.loc[metadata['passage']>0,:].sort_values(by='passage')
    n_snps = len(parent_snps)

    diff_info_lower=[]
    diff_info_upper=[]
    diff_info_med = []
    pred_7_upper = []
    pred_7_lower= []
    pred_med = []
    mesocosms = []
    for mesocosm in metadata['mesocosm'].unique():
        meso_df = metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'sample'].values
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values)
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        if inoculumn not in metadata['sample'].values:
            continue 
        if 3 not in passages:
            continue
        if 7 not in passages:
            continue 
            
        sample0 = inoculumn 
        sample3 = meso_df.loc[meso_df['passage']==3,'sample'].values
        sample7 = meso_df.loc[meso_df['passage']==7,'sample'].values
        
        diffs = np.zeros(n_bootstraps)
        pred_7_ests = np.zeros(n_bootstraps)
        est_7s_ests = np.zeros(n_bootstraps)
        for sample in range(n_bootstraps):
            bs_sample = np.random.shuffle(parent_snps)
            int1 = round(len(bs_sample)/3)
            int3 = round(2*len(bs_sample)/3)

            passage0_est = get_freq_est(snps.loc[bs_sample[:int1],sample0].values, 
                        depth_med[sample0] )
            passage3_est = get_freq_est(snps.loc[bs_sample[int1+1:int2],sample3].values, 
                        depth_med[sample0] )
            passage7_est = get_freq_est(snps.loc[bs_sample[int2+1:],sample7].values, 
                        depth_med[sample0] )
            
            dt = 3
            sel_coeff_03 = (1/dt)*(np.log(passage3_est/(1-passage3_est)) + \
                                np.log((1-passage0_est)/(passage0_est)))
            pred_p7 = np.exp(7*sel_coeff_03)
            pred_p7 = pred_p7*passage0_est/(1 - passage0_est + passage0_est*pred_p7)

            diffs[sample] = pred_p7 - passage7_est 
            est_7s_ests[sample] = passage7_est
            pred_7_ests[sample] = pred_p7 
        diff_info_upper.append(np.quantile(diffs q = .975))
        diff_info_lower.append(np.quantile(diffs q = .025))
        pred_7_upper.append(np.quantile(pred_7_ests, q = .975))
        pred_7_lower.append(np.quantile(pred_7_ests, q = .025))
        pred_med.append(np.median(pred_7_ests))
        diff_info_med.append(np.median(diffs))
        mesocosms.append(mesocosm)
    df = pd.DataFrame(data = { 'mesocosms': mesocosms, 'diff_med': diff_info_med, 'diff_lower': diff_info_lower, 
                                'diff_upper': diff_info_upper, 'pred_7_upper': pred_7_upper,
                                'pred_7_lower': pred_7_lower, 'pred_med': pred_med,
                                })
    df['n_snps'] = n_snps
    return df
       



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

def get_main(metadata, species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered):
  #  info, depth, freq = load_and_sort_files(species_dir, species)
   # med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
   # good_samples = med_nonzero_depth[med_nonzero_depth>10.]
   # depth = depth[good_samples.index.values]
   # freq = freq[good_samples.index.values]  
   # depth_filtered= depth_filtering(depth)
   # freq_filtered = freq_masked(freq, depth_filtered)

    freq_inoculumns = freq_filtered[parent_samples]
    print(parent_samples)
    
    if len(freq_inoculumns.columns.values)<2:
        return pd.DataFrame(), pd.DataFrame()
    # get distinguishing SNPs for inoculumns - there should be like 1k distinguishing SNPs 
    # use only Alt Allele as marker... so sites where strain allele is alt allele in one strain and not other strain
    # is the marker 
    distinguishing_snps = get_distinguishing_snps(freq_inoculumns, thresh = .8)
   # print(len(distinguishing_snps))
    #distinguishing_snps.to_csv('distinguishing_snps.csv')
    parent1_snps = distinguishing_snps[distinguishing_snps == 1].index.values
    parent2_snps = distinguishing_snps[distinguishing_snps == -1].index.values
    freq_children = freq_filtered[child_samples]
    depth_children = depth_filtered[child_samples]
   
    print('before', len(parent1_snps), len(parent2_snps))
    parent1_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent1_snps, thresh = .5, sample_thresh=.75)
    parent2_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent2_snps, thresh = .5, sample_thresh=.75)
    print('after', len(parent1_snps), len(parent2_snps))
    med_depth_children = depth_children.copy().replace(0, np.nan).median(skipna=True)

   # count_parent1 = pd.DataFrame(get_count(freq_children, parent1_snps)).rename(columns = {0: parent_samples[0]})  
   # count_parent2 = pd.DataFrame(get_count(freq_children, parent2_snps)).rename(columns = {0: parent_samples[1]})  

    parent1_info = get_bootstrap_sel_coeffs(metadata, freq_children, med_depth_children, parent1_snps, n_bootstraps = 1000)
    parent2_info = get_bootstrap_sel_coeffs(metadata, freq_children, med_depth_children, parent2_snps, n_bootstraps = 1000)
# freq_parent2 = freq_parent2.T
   # freq_parent2['parent'] = parent_samples[1]
   # print(freq_parent2)
    return parent1_info, parent2_info 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
#    parser.add_argument('--inoculumn', action = 'store', 
 #                      help = 'inoculumn')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'

#    parent_samples, child_samples = get_parent_children(args.inoculumn)


    if not path.isdir(save_dir):
        mkdir(save_dir) 
    info, depth, freq = load_and_sort_files(species_dir, args.species)
    #print(info.columns.values)
    #print(info.index.values)
    freq = repolarize_against_reference(freq, info)
    metadata = pd.read_csv('workflow/analysis/e003_metadata_cultures_round2_change_AA.csv')

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    for inoculumn in inoculumn_list:
        print(inoculumn)
        parent_samples, child_samples = get_parent_children(inoculumn, metadata)
      #  print(parent_samples)
#        print(freq_filtered.columns.values)
#AAAA
        parent_samples = list(np.intersect1d(parent_samples, freq_filtered_in.columns.values))
        child_samples = list(np.intersect1d(child_samples, freq_filtered_in.columns.values))
        #child_samples = freq_filtered_in.columns.values
        print(parent_samples)
        print(child_samples)
        if len(parent_samples) < 2:
            continue
        if parent_samples[1] in child_samples:
            child_samples.remove(parent_samples[1])
        if parent_samples[0] in child_samples:
            child_samples.remove(parent_samples[0])
        if parent_samples[0] == 'A2-e003Coalescence-Inoculumn-mBHI':
            parent_samples = ['A2-e003Coalescence-mBHI-inoculumn-redo', parent_samples[1]]
        
        parent1_info, parent2_info = get_main(metadata,species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples])

        inoculumn = ''.join(inoculumn.split('/'))
    
      #  count_parents.to_csv(f'{save_dir}/{inoculumn}_count_parents.csv')
        parent1_info.to_csv(f'{save_dir}/{inoculumn}_parent1_info.csv')
        parent2_info.to_csv(f'{save_dir}/{inoculumn}_parent2_info.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

