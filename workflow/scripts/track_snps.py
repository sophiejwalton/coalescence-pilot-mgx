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
#    print(detect_df[freq_inoculumns.columns.values[0]])
    detect_df['diff'] = detect_df[freq_inoculumns.columns.values[0]].astype(int) - detect_df[freq_inoculumns.columns.values[1]].astype(int)
    return detect_df['diff']

def get_frequency_parent(freq_children, parent_snps):
    median_freq = freq_children.loc[parent_snps,:].median(axis = 0)

    return median_freq 


def get_parent_children(inoculumn):
    metadata = pd.read_csv('config/e003_coal_metadata_full.csv')
    child_samples = list(metadata.loc[metadata['inoculumn'] == inoculumn, 'sample'].values)
    parent_subjects = inoculumn.split('-')[:-1]
    parent_media = inoculumn.split('-')[-1]
    ins = metadata.loc[metadata['is_inoculumn'],:]
    ins = ins.loc[ins['parent_media'] == parent_media,:]
    ins1 = ins.loc[ins['parent_subjects'] == parent_subjects[0] + '-' +  parent_subjects[0],:]

    ins2 = ins.loc[ins['parent_subjects'] == parent_subjects[1] + '-' +  parent_subjects[1],:]
    if (len(ins1) ==0) or (len(ins2) ==0): 
        return np.nan, np.nan
    parent_samples = [ins1['sample'].values[0], ins2['sample'].values[0]]
    return parent_samples, child_samples
    


def get_main(species_dir,species, parent_samples, child_samples, freq_filtered):
  #  info, depth, freq = load_and_sort_files(species_dir, species)
   # med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
   # good_samples = med_nonzero_depth[med_nonzero_depth>10.]
   # depth = depth[good_samples.index.values]
   # freq = freq[good_samples.index.values]  
   # depth_filtered= depth_filtering(depth)
   # freq_filtered = freq_masked(freq, depth_filtered)

    freq_inoculumns = freq_filtered[parent_samples]
    
    if len(freq_inoculumns.columns.values)<2:
        return pd.DataFrame()
    # get distinguishing SNPs for inoculumns - there should be like 1k distinguishing SNPs 
    # use only Alt Allele as marker... so sites where strain allele is alt allele in one strain and not other strain
    # is the marker 
    distinguishing_snps = get_distinguishing_snps(freq_inoculumns, thresh = .999)
    print(distinguishing_snps)
    distinguishing_snps.to_csv('distinguishing_snps.csv')
    parent1_snps = distinguishing_snps[distinguishing_snps == 1].index.values
    parent2_snps = distinguishing_snps[distinguishing_snps == -1].index.values
   # print(parent1_snps)
    #print(parent2_snps)
    freq_children = freq_filtered[child_samples]

    freq_parent1 = get_frequency_parent(freq_children, parent1_snps)
#    print(freq_parent1)
    freq_parent1 = pd.DataFrame(freq_parent1).rename(columns = {0: parent_samples[0]})
   # print(freq_parent1)
    freq_parent1.to_csv('freq_parent1.csv')
  #  freq_parent1['parent'] = parent_samples[0]
    freq_parent2 = get_frequency_parent(freq_children, parent2_snps)
    freq_parent2 = pd.DataFrame(freq_parent2).rename(columns = {0: parent_samples[1]})  
# freq_parent2 = freq_parent2.T
   # freq_parent2['parent'] = parent_samples[1]
   # print(freq_parent2)
    return pd.concat([freq_parent1, freq_parent2],axis=1) 

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
    print(info.columns.values)
    print(info.index.values)
    freq = repolarize_against_reference(freq, info)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>10.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    for inoculumn in inoculumn_list:
        print(inoculumn)
        parent_samples, child_samples = get_parent_children(inoculumn)
        parent_samples = list(np.intersect1d(parent_samples, freq_filtered.columns.values))
        child_samples = list(np.intersect1d(child_samples, freq_filtered.columns.values))
        _, freq_filtered_in = filter_sites_across_samples(depth_filtered[parent_samples + child_samples], freq_filtered[parent_samples + child_samples].copy(), )

        freq_parents = get_main(species_dir,args.species, parent_samples, child_samples, freq_filtered_in)
        inoculumn = ''.join(inoculumn.split('/'))
        freq_parents.to_csv(f'{save_dir}/{inoculumn}_parent_freqs.csv')
    


       
        
    

