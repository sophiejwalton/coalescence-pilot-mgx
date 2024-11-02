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


def get_haplotype(sample_freq_filtered,sample):
    good_sites_lower = sample_freq_filtered.loc[sample_freq_filtered[sample] < .2]
    good_sites_upper = sample_freq_filtered.loc[sample_freq_filtered[sample] > .8]
    all_good_sites = pd.concat([good_sites_lower, good_sites_upper])
    all_good_sites = all_good_sites.loc[~np.isnan(all_good_sites[sample]),:]
    all_good_sites['Allele'] = 1*(all_good_sites[sample] > .5)
    return all_good_sites




def get_pairwise_divergence(sample1_haplotype, sample2_haplotype):
   # print(sample1_haplotype, sample2_haplotype) 
    shared_sites = np.intersect1d(sample1_haplotype.index.values, sample2_haplotype.index.values)
    sample1_haplotype_shared_site = sample1_haplotype.loc[shared_sites, :].copy()
    sample2_haplotype_shared_site = sample2_haplotype.loc[shared_sites, :].copy()
    same_allele = sample1_haplotype_shared_site['Allele'] == sample2_haplotype_shared_site['Allele']
    
    fixed_differences = len(same_allele) - same_allele.sum()
    divergence = fixed_differences/len(same_allele)
    return fixed_differences, divergence, len(same_allele)
    


good_inoculumns = {'AA-AA-mBHI': ['A2-e003Coalescence-mBHI-inoculumn-redo',
  'A2-e003Coalescence-Inoculumn-mBHI'],
 'AA-AA-mGAM': ['A2-e003Coalescence-Inoculumn-mGAM'],
 'AA-AC/PP-mBHI': ['A3-e003Coalescence-mBHI-inoculumn-redo'],
 'AA-AC/PP-mGAM': ['A3-e003Coalescence-Inoculumn-mGAM'],
 'AA-AE-mBHI': ['A4-e003Coalescence-mBHI-inoculumn-redo'],
 'AA-AE-mGAM': ['A4-e003Coalescence-Inoculumn-mGAM'],
 'AA-AF-mBHI': ['A5-e003Coalescence-mBHI-inoculumn-redo',
  'A5-e003Coalescence-Inoculumn-mBHI'],
 'AA-AF-mGAM': ['A5-e003Coalescence-Inoculumn-mGAM'],
 'AC/PP-AA-mBHI': ['B2-e003Coalescence-Inoculumn-mBHI'],
 'AC/PP-AA-mGAM': ['B2-e003Coalescence-Inoculumn-mGAM'],
 'AC/PP-AC/PP-mBHI': ['B3-e003Coalescence-Inoculumn-mBHI'],
 'AC/PP-AC/PP-mGAM': ['B3-e003Coalescence-Inoculumn-mGAM'],
 'AC/PP-AE-mBHI': ['B4-e003Coalescence-Inoculumn-mBHI'],
 'AC/PP-AE-mGAM': ['B4-e003Coalescence-Inoculumn-mGAM'],
 'AC/PP-AF-mBHI': ['B5-e003Coalescence-Inoculumn-mBHI'],
 'AC/PP-AF-mGAM': ['B5-e003Coalescence-Inoculumn-mGAM'],
 'AE-AA-mBHI': ['C2-e003Coalescence-mBHI-inoculumn-redo'],
 'AE-AA-mGAM': ['C2-e003Coalescence-Inoculumn-mGAM'],
 'AE-AC/PP-mBHI': ['C3-e003Coalescence-mBHI-inoculumn-redo'],
 'AE-AC/PP-mGAM': ['C3-e003Coalescence-Inoculumn-mGAM'],
 'AE-AE-mBHI': ['C4-e003Coalescence-mBHI-inoculumn-redo'],
 'AE-AE-mGAM': ['C4-e003Coalescence-Inoculumn-mGAM'],
 'AE-AF-mBHI': ['C5-e003Coalescence-mBHI-inoculumn-redo'],
 'AE-AF-mGAM': ['C5-e003Coalescence-Inoculumn-mGAM'],
 'AF-AA-mBHI': ['D2-e003Coalescence-Inoculumn-mBHI'],
 'AF-AA-mGAM': ['D2-e003Coalescence-Inoculumn-mGAM'],
 'AF-AC/PP-mBHI': ['D3-e003Coalescence-Inoculumn-mBHI'],
 'AF-AC/PP-mGAM': ['D3-e003Coalescence-Inoculumn-mGAM'],
 'AF-AE-mBHI': ['D4-e003Coalescence-Inoculumn-mBHI'],
 'AF-AE-mGAM': ['D4-e003Coalescence-Inoculumn-mGAM'],
 'AF-AF-mBHI': ['D5-e003Coalescence-Inoculumn-mBHI'],
 'AF-AF-mGAM': ['D5-e003Coalescence-Inoculumn-mGAM']}

 
def get_main(species_dir, save_dir, species, parent_subjects_media):
    parent_subjects_media = parent_subjects_media.split('-')
    media = parent_subjects_media[-1]
    parent_subject = parent_subjects_media[:-1]
    info, depth, freq = load_and_sort_files(species_dir, species)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]  
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)

    inoculumns = good_inoculumns[f'{parent_subject[0]}-{media}']
 #   freq_inoculumns = freq_filtered[inoculumns]

   # = get_haplotype(freq_filtered, inoculumns[0])

   # haplotype2 = get_haplotype(freq_filtered, inoculumns[1])


    for i, s1 in enumerate(good_inoculumns):
        print(i, s1)
        if s1 not in good_samples.index.values:
            continue 

        for s2 in good_samples.index.values: 
            print(s2)
            if s1 == s2:
                continue 
            freq_small = freq_filtered[[s1,s2]].copy()
            depth_small =  depth_filtered[[s1,s2]].copy()
            freq_small = polarize_species(freq_small.copy(), s1)
            freq_polarized_transition = get_transition_frequency_snps(freq_small, depth_small)
            depth_small = depth_small.loc[freq_polarized_transition.index.values, ]
            freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(), depth_small, med_nonzero_depth[[s1,s2]])
            snps_switch.append(len(freq_transition_filter))
            s1s.append(s1)
            s2s.append(s2)
    ss_df = pd.DataFrame(data = {'s1': s1s, 's2': s2s, 'fixed_diffs': snps_switch, })
    ss_df['Species'] = species
    ss_df['Strain Shift'] = ss_df['fixed_diffs'] > 1000
    if '/' in parent_subjects_media:
        parent_subjects_media = ''.join(parent_subjects_media.split('/'))
    ss_df.to_csv(f'{save_dir}/{species}_{parent_subjects_media}_fixed_diffs.csv')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
    parser.add_argument('--parent_subject-media', action = 'store', 
                       help = 'parent_subject-media')
    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species, args.parent_subject-media)
    


       
        
    

