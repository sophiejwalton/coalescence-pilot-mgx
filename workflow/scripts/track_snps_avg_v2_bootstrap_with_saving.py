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



def get_main(species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered, samples_to_save):
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

    freq_children.loc[parent1_snps, samples_to_save].to_csv(f'{species_dir}/to_save_dist1.csv.gz',compression='gzip')
    freq_children.loc[parent2_snps, samples_to_save].to_csv(f'{species_dir}/to_save_dist2.csv.gz',compression='gzip')


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
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    #freq_filtered_in[['D_B8_e003Assembly_S308' ]].to_csv(f'{species_dir}/for_plotting_histogram.csv.gz',compression='gzip')
    #metadata = pd.read_csv('workflow/analysis/e003_metadata_cultures_round2_change_AA.csv')
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4.csv')
    for inoculumn in ['AE-AF-mBHI']:
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

        
        samples_to_save = ['C5-e003Coalescence-mBHI-inoculumn-redo','B_G4_G5_AF_AE_mBHI_mGAM_2_S172',
       'C_G12_G5_AF_AE_mBHI_mGAM_4_S276',
       'D_G4_C5_AF_AE_mBHI_mGAM_6_S364',
       'G5-e003Coalescence-mGAM-p1_S810',
       'G5-e003Coalescence-mGAM-p3_S84', 'G5-e003Coalescence-mGAM-p5_S76',
       'G5-e003Coalescence-mGAM-p7_S172']
        samples_to_save = ['A11-AA-AF-mGAM-mBHI', 'A11-AA-AF-mGAM-mGAM', 'A5-AA-AF-mGAM-mGAM',
 'D8-AA-AF-mGAM-mBHI', 'D8-AA-AF-mGAM-mGAM', 'E11-AA-AF-mGAM-mBHI',
 'H8-AA-AF-mGAM-mBHI', 'E11-AA-AF-mGAM-mGAM', 'H8-AA-AF-mGAM-mGAM',
 'D2-AA-AF-mGAM-mGAM']
        samples_to_save=['Coalescence-F1-D8-AA-AF-mGAM-mGAM-3_S680',
       'D8-e003Coalescence-mGAM-p5',
       'Coalescence-B5-D8-AA-AF-mGAM-mGAM-7_S636',
       'C_H5_D8_AA_AF_mGAM_mGAM_6_S281', 'A_H5_D8_AA_AF_mGAM_mGAM_2_S89',
       'B_H9_D8_AA_AF_mGAM_mGAM_4_S189',
       'D8-e003Coalescence-mGAM-p1_S755',
       'A5-e003Coalescence-Inoculumn-mGAM']
        samples_to_save = ['B_E8_E11_AF_AA_mGAM_mGAM_2_S152',
       'Coalescence-B6-E11-AA-AF-mGAM-mGAM-7_S637',
        'A5-e003Coalescence-Inoculumn-mGAM',
       'Coalescence-F3-E11-AA-AF-mGAM-mGAM-3_S682',
       'D_E12_E11_AF_AA_mGAM_mGAM_4_S348',
       'D_E8_A11_AF_AA_mGAM_mGAM_6_S344',
       'E11-e003Coalescence-mGAM-p1_S840',
       'E111-e003Coalescencei-mGAM-p5_S56']
        samples_to_save = ['A_C12_C5_AF_AE_mBHI_mBHI_4_S36', 'A_C4_C5_AF_AE_mBHI_mBHI_2_S28',
       'C5-e003Coalescence-mBHI-inoculumn-redo',
       'C5-e003Coalescence-mBHI-p3', 'C5-e003Coalescence-mBHI-p5',
       'C5-e003Coalescence-mBHI-p7', 'C_C4_C5_AF_AE_mBHI_mBHI_6_S220']
        samples_to_save = ['D_C4_C5_AF_AE_mBHI_mBHI_6_S316',
       'G5-e003Coalescence-mBHI-p7_S124',
       'C_C12_G5_AF_AE_mBHI_mBHI_4_S228',
       'G5-e003Coalescence-mBHI-p5_S28', 'B_C4_G5_AF_AE_mBHI_mBHI_2_S124',
       'G5-e003Coalescence-mBHI-p3_S36',
       'G5-e003Coalescence-mBHI-p1_S806',
       'C5-e003Coalescence-mBHI-inoculumn-redo']
        get_main(species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples],samples_to_save)

       
        
    

