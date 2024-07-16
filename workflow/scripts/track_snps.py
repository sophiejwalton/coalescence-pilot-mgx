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
    detect_df = freq_inoculumns > thresh
    print(detect_df[freq_inoculumns.columns.values[0]])
    detect_df['diff'] = detect_df[freq_inoculumns.columns.values[0]] - detect_df[freq_inoculumns.columns.values[1]]
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
    


def get_main(species_dir,species, parent_samples, child_samples):
    info, depth, freq = load_and_sort_files(species_dir, species)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>10.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)

    freq_inoculumns = freq_filtered[parent_samples]
    # get distinguishing SNPs for inoculumns - there should be like 1k distinguishing SNPs 
    # use only Alt Allele as marker... so sites where strain allele is alt allele in one strain and not other strain
    # is the marker 
    distinguishing_snps = get_distinguishing_snps(freq_inoculumns, thresh = .8)
    parent1_snps = distinguishing_snps.loc[distinguishing_snps == 1,:].index.values
    parent2_snps = distinguishing_snps.loc[distinguishing_snps == -1,:].index.values

    freq_children = freq_filtered[child_samples]

    freq_parent1 = get_frequency_parent(freq_children, parent1_snps)
    freq_parent1 = freq_parent1.T
    freq_parent1['parent'] = parent_samples[0]
    freq_parent2 = get_frequency_parent(freq_children, parent2_snps)
    freq_parent1 = freq_parent1.T
    freq_parent1['parent'] = parent_samples[1]
    return pd.concat([freq_parent1, freq_parent2])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
    parser.add_argument('--inoculumn', action = 'store', 
                       help = 'inoculumn')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'

    parent_samples, child_samples = get_parent_children(args.inoculumn)


    if not path.isdir(save_dir):
        mkdir(save_dir)


    freq_parents = get_main(species_dir,args.species, parent_samples, child_samples)
    freq_parents.to_csv(f'{save_dir}/{inoculumn}_parent_freqs.csv')
    


       
        
    

