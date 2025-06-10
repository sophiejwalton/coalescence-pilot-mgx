import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
import itertools as it
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
from track_snps_funcs import * 
import warnings
warnings.filterwarnings('ignore')


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
    passage1s=[]
    passage2s=[]
    sample1s=[]
    sample2s=[]
    sel_inf_lower=[]
    sel_inf_upper=[]
    sel_inf_med = []
    mesocosms = []
    for mesocosm in metadata['mesocosm'].unique():
        samples = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'sample'].values)
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values)
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        if inoculumn not in metadata['sample'].values:
            continue 
        samples = [inoculumn] + samples
        passages = [0] + passages
        print(samples, snps.columns.values)
       #snps_samples = snps[samples].values 
       # depths = depth_med[samples].values
        medians = []
        for sample in samples:
            snps_sample = snps[sample].values
            snps_sample = snps_sample[~np.isnan(snps_sample)]
            n_snps = len(snps_sample)
            bs_samples = np.random.choice(snps_sample, size = ( len(snps), n_bootstraps))
            meds = np.nanmedian(bs_samples, axis=0)

            b_frac_zero = np.sum(bs_samples == 0, axis = 0)/n_snps
            b_frac_zero = -np.log(b_frac_zero)/depth_med[sample]

            b_frac_one = np.sum(bs_samples == 1, axis = 0)/n_snps
            b_frac_one = 1+ np.log(b_frac_one)/depth_med[sample]

            meds[meds == 0 ] = b_frac_zero[meds == 0]
            meds[meds == 1 ] = b_frac_one[meds == 1]
            medians.append(meds)

        for combo in it.combinations(np.arange(len(passages)),2):
            combo = np.sort([combo[0], combo[1]])
            passage1 = passages[combo[0]]
            passage2 = passages[combo[1]]
            sample1 = samples[combo[0]]
            sample2 = samples[combo[1]]
            dt = passage2-passage1

          #  all_medians = medians[:,combo]
            sel_coeffs = (1/dt)*(np.log(medians[combo[1]]/(1-medians[combo[1]])) + \
                             np.log((1-medians[combo[0]])/(medians[combo[0]])))
            sel_inf_upper.append(np.quantile(sel_coeffs, q = .975))
            sel_inf_lower.append(np.quantile(sel_coeffs,q=.025))
            sel_inf_med.append(np.median(sel_coeffs))
            passage1s.append(passage1)
            passage2s.append(passage2)
            sample1s.append(sample1)
            sample2s.append(sample2)
            mesocosms.append(mesocosm)
    df = pd.DataFrame(data = {'sample1': sample1s, 'sample2': sample2s, 'passage1': passage1s, 'passage2': passage2s, 
                                'mesocosms': mesocosms, 'sel_med': sel_inf_med, 'sel_lower': sel_inf_lower, 
                                'sel_upper': sel_inf_upper, 
                                })
    df['n_snps'] = n_snps
    return df
       


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
    distinguishing_snps1, distinguishing_snps2= get_distinguishing_snps(freq_inoculumns, thresh = .8)
   # print(len(distinguishing_snps))
    #distinguishing_snps.to_csv('distinguishing_snps.csv')
    parent1_snps = distinguishing_snps1[distinguishing_snps1 == 2].index.values
    parent2_snps = distinguishing_snps2[distinguishing_snps2 == 2].index.values
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
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4.csv')

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
        parent_samples, child_samples = get_parent_children(inoculumn,metadata)
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
    


       
        
    

