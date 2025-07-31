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

def adjust_freq(nreads0, nreads3):
    return (6*nreads3/nreads0)**(1/3)

def get_sample_freq(freqs, reads,readsopp, depth_med):
    reads = reads[~np.isnan(freqs)]
    readsopp = readsopp[~np.isnan(freqs)]
    freqs = freqs[~np.isnan(freqs)]
    freq_est = np.median(freqs)
    if freq_est == 0:
        freq_est = adjust_freq(np.sum(reads==0), np.sum(reads==3))/depth_med
     #   if np.sum(reads == 3) > 100:
      #      if np.sum(reads == 4)/np.sum(reads == 3) <.9
       #         freq_est = (np.sum(reads==4)*4/np.sum(reads==3))/depth_med
    if freq_est == 1:
        freq_est = 1-adjust_freq(np.sum(readsopp==0), np.sum(readsopp==3))/depth_med
       # if np.sum(readsopp == 3) > 100:
        #    if np.sum(readsopp == 4)/np.sum(readsopp == 3) <.9
         #       freq_est = (np.sum(readsopp==4)*4/np.sum(readsopp==3))/depth_med
    return freq_est 

def get_bootstrap_parent(metadata,freq_children, depth_med, depth_children, reads, readsopp, parent_snps1,parent_snps2, inoculumn,thresh = 1e-2, n_bootstraps = 1000):
    passage1s=[]
    passage2s=[]
    sample1s=[]
    sample2s=[]
    sel_inf_lower=[]
    sel_inf_upper=[]
    sel_inf_med = []
    mesocosms = []
    shifts1=[]
    shifts2=[]
    initial_freqs=[]
    metadata = metadata.loc[metadata['sample'].isin(freq_children.columns.values),:]
    metadata_non_zero = metadata.loc[metadata['passage']>0,:].sort_values(by='passage')
    for mesocosm in metadata['mesocosm'].unique():
        samples = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'sample'].values)
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values)
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        if inoculumn not in metadata['sample'].values:
            continue 
        samples = [inoculumn] + samples
        passages = [0] + passages
        medians = []
        for combo in it.combinations(np.arange(len(passages)),2):
            combo = np.sort([combo[0], combo[1]])
            passage1 = passages[combo[0]]
            passage2 = passages[combo[1]]
            sample1 = samples[combo[0]]
            sample2 = samples[combo[1]]
            dt = passage2-passage1
            sel_coeffs = np.zeros(n_bootstraps)
            shifts_pair1 = np.zeros(n_bootstraps)
            shifts_pair2 = np.zeros(n_bootstraps)
            initial_freqs_pair = np.zeros(n_bootstraps)
            for bs in range(n_bootstraps):
                bs_strain1 = np.random.choice(parent_snps1, size=len(parent_snps1),replace=False)
                bs_strain2 = np.random.choice(parent_snps2, size=len(parent_snps2),replace=False)
                snps_strain1_sample1 = bs_strain1[:int(len(parent_snps1)/2)]
                snps_strain2_sample1 = bs_strain2[:int(len(parent_snps2)/2)]
                freq_strain1_sample1 = get_sample_freq(freq_children.loc[snps_strain1_sample1,:], reads.loc[snps_strain1_sample1,:],readsopp.loc[snps_strain1_sample1,:], depth_med[sample1])
                freq_strain2_sample1 = get_sample_freq(freq_children.loc[snps_strain2_sample1,:], reads.loc[snps_strain2_sample1,:],readsopp.loc[snps_strain2_sample1,:], depth_med[sample1])

                snps_strain1_sample2 = bs_strain1[int(len(parent_snps1)/2)+1:]
                snps_strain2_sample2 = bs_strain2[int(len(parent_snps2)/2)+1:]
                freq_strain1_sample2 = get_sample_freq(freq_children.loc[snps_strain1_sample2,:], reads.loc[snps_strain1_sample2,:],readsopp.loc[snps_strain1_sample2,:], depth_med[sample2])
                freq_strain2_sample2 = get_sample_freq(freq_children.loc[snps_strain2_sample2,:], reads.loc[snps_strain2_sample2,:],readsopp.loc[snps_strain2_sample2,:], depth_med[sample2])

                shifts_pair1[bs] = np.abs(freq_strain1_sample1-freq_strain2_sample1)
                shifts_pair2[bs] = np.abs(freq_strain1_sample2-freq_strain2_sample2)
                freq_sample1 = freq_strain1_sample1/(freq_strain1_sample1+freq_strain2_sample1 )
                if freq_sample1<thresh:
                    freq_sample1 = thresh
                if freq_sample1>1-thresh:
                    freq_sample1 = 1-thresh

                freq_sample2 = freq_strain1_sample2/(freq_strain1_sample2+freq_strain2_sample2 )
                initial_freqs_pair[bs]=freq_sample1

                if freq_sample2<thresh:
                    freq_sample2 = thresh
                if freq_sample2>1-thresh:
                    freq_sample2 = 1-thresh
                sel_coeffs[bs]= (1/dt)*(np.log(freq_sample2/(1-freq_sample2)) + \
                             np.log((1-freq_sample1)/(freq_sample1)))
            passage1s.append(passage1)
            passage2s.append(passage2)
            sample1s.append(sample1)
            sample2s.append(sample1)
            mesocosms.append(mesocosm)
            sel_inf_lower.append(np.percentile(sel_coeffs, q =2.5))
            sel_inf_upper.append(np.percentile(sel_coeffs,q=97.5))
            sel_inf_med.append(np.median(sel_coeffs))
            shifts1.append(np.median(shifts_pair1))
            shifts2.append(np.median(shifts_pair2))
            initial_freqs.append(np.median(initial_freqs_pair))
            
    return  pd.DataFrame(data = {'sample1': sample1s, 'sample2': sample2s, 'passage1': passage1s, 'passage2': passage2s, 
                                'mesocosms': mesocosms, 'sel_med': sel_inf_med, 'sel_lower': sel_inf_lower, 
                                'sel_upper': sel_inf_upper, 'shifts1': shifts1, 'shifts2': shifts2, 'initial_freqs': initial_freqs
                                })


def get_main(metadata,species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered,inoculumn):
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
    reads_children = freq_children*depth_children
    reads_children_opp = (1-freq_children)*depth_children

   # print(reads_children)
   
    print('before', len(parent1_snps), len(parent2_snps))
    parent1_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent1_snps, thresh = .5, sample_thresh=.75)
    parent2_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent2_snps, thresh = .5, sample_thresh=.75)
    print('after', len(parent1_snps), len(parent2_snps))
    med_depth_children = depth_children.copy().replace(0, np.nan).median(skipna=True)

    count_parent1 = pd.DataFrame(get_count(freq_children, parent1_snps)).rename(columns = {0: parent_samples[0]})  
    count_parent2 = pd.DataFrame(get_count(freq_children, parent2_snps)).rename(columns = {0: parent_samples[1]})  

    parent1_info = get_bootstrap_parent(metadata, freq_children, med_depth_children, depth_children, reads_children, reads_children_opp, parent1_snps,parent2_snps, inoculumn, n_bootstraps = 1000)
    #parent2_info = get_bootstrap_parent(freq_children, med_depth_children,depth_children, reads_children, reads_children_opp, parent2_snps, inoculumn, n_bootstraps = 1000)
#    reads_children.loc[reads_children,parent1_snps].to_csv(f'{species_dir}/{inoculumn}_reads_children.csv.gz',compression='gzip')
 #   reads_children_opp.loc[reads_children,parent1_snps].to_csv(f'{species_dir}/{inoculumn}_reads_children_opp.csv.gz',compression='gzip')
    return pd.concat([count_parent1, count_parent2],axis=1), parent1_info

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

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
        'AA-AE-mBHI', 'AA-AF-mBHI',
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    #metadata = pd.read_csv('workflow/analysis/e003_metadata_cultures_round2_change_AA.csv')
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4_good.csv')
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
        
        count_parents, parent1_info = get_main(metadata,species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples],inoculumn)

        inoculumn = ''.join(inoculumn.split('/'))
        

        count_parents.to_csv(f'{save_dir}/{inoculumn}_count_parents.csv')
        parent1_info.to_csv(f'{save_dir}/{inoculumn}_parent1_info.csv')
       # parent2_info.to_csv(f'{save_dir}/{inoculumn}_parent2_info.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

