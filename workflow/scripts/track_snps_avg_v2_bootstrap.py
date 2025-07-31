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
        if np.sum(reads == 3) > 100:
            if np.sum(reads == 4)/np.sum(reads == 3) <.9:
                freq_est = (np.sum(reads==4)*4/np.sum(reads==3))/depth_med
    if freq_est == 1:
        freq_est = 1-adjust_freq(np.sum(readsopp==0), np.sum(readsopp==3))/depth_med
        if np.sum(readsopp == 3) > 100:
            if np.sum(readsopp == 4)/np.sum(readsopp == 3) <.9:
                freq_est = (np.sum(readsopp==4)*4/np.sum(readsopp==3))/depth_med
    return freq_est 

def get_bootstrap_parent(freq_children, depth_med,depth_children, reads_children, reads_children_opp, parent_snps, inoculumn,n_bootstraps = 1000):
    snps = freq_children.loc[parent_snps,:]
    depths = depth_children.loc[parent_snps,:]
    reads_opp=reads_children_opp.loc[parent_snps,:]
    boot_med = []
    boot_low = []
    boot_high = []
    act_med = []
    samples_good = []
    for sample in snps.columns.values:
        snps_sample= snps[sample]
        reads_sample = depths[sample]*snps_sample
 #       print(reads_sample,snps_sample)
        reads_sample = reads_sample[~np.isnan(snps_sample)]
  #      print(np.min(reads_sample[reads_sample>0]))
        reads_sample = reads_sample.round()
   #     print(np.min(reads_sample[reads_sample>0]))
        reads_sample_opp = (depths[sample]*(1-snps_sample))
        reads_sample_opp = reads_sample_opp[~np.isnan(snps_sample)].round()
        snps_sample = snps_sample[~np.isnan(snps_sample)]
        med_og = get_sample_freq(snps_sample, reads_sample,reads_sample_opp, depth_med[sample])
        if len(snps_sample) == 0: continue 
     #   print('after',med_og)
#        bs_samples = np.random.choice(snps_sample.index.values, size = (n_snps, n_bootstraps))
        samples_meds= np.zeros(n_bootstraps)
        for n in range(n_bootstraps):
            bs_sample = np.random.choice(snps_sample.index.values, size = len(snps_sample))
            snps_bs = snps_sample[bs_sample]
            reads_bs = reads_sample[bs_sample]
            reads_opp_bs = reads_sample_opp[bs_sample]
            bs_med = get_sample_freq(snps_bs, reads_bs,reads_opp_bs, depth_med[sample])
            samples_meds[n] = bs_med
            
        samples_good.append(sample)
        boot_med.append(np.median(samples_meds))
        boot_low.append(np.percentile(samples_meds, q =2.5))
        act_med.append(med_og)
        boot_high.append(np.percentile(samples_meds, q = 97.5))
    return pd.DataFrame(data = {'sample': samples_good, 'boot_med': boot_med, 'boot_low': boot_low, 'boot_high': boot_high,
                                    'actual_med': act_med})


def get_main(species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered,inoculumn):
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

    parent1_info = get_bootstrap_parent(freq_children, med_depth_children,depth_children, reads_children, reads_children_opp, parent1_snps,inoculumn, n_bootstraps = 1000)
    parent2_info = get_bootstrap_parent(freq_children, med_depth_children,depth_children, reads_children, reads_children_opp, parent2_snps, inoculumn, n_bootstraps = 1000)
#    reads_children.loc[reads_children,parent1_snps].to_csv(f'{species_dir}/{inoculumn}_reads_children.csv.gz',compression='gzip')
 #   reads_children_opp.loc[reads_children,parent1_snps].to_csv(f'{species_dir}/{inoculumn}_reads_children_opp.csv.gz',compression='gzip')
    return pd.concat([count_parent1, count_parent2],axis=1), parent1_info, parent2_info 

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
        
        count_parents, parent1_info, parent2_info = get_main(species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples],inoculumn)

        inoculumn = ''.join(inoculumn.split('/'))
        

        count_parents.to_csv(f'{save_dir}/{inoculumn}_count_parents.csv')
        parent1_info.to_csv(f'{save_dir}/{inoculumn}_parent1_info.csv')
        parent2_info.to_csv(f'{save_dir}/{inoculumn}_parent2_info.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

