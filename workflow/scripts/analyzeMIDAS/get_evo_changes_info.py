import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
from tqdm import tqdm
import argparse

from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
warnings.filterwarnings('ignore')


def get_main(subjects, species_dir, save_dir, species,  main_study = True):
    info, _, _ = load_and_sort_files(species_dir)
    info = info.set_index('site_id')
    num_genes_hit = np.zeros(len(subjects))
    num_snps_change=  np.zeros(len(subjects))
    for i, subject in enumerate(subjects):
        if args.main_study:
            freq_transition = pd.read_csv(f'{save_dir}/{subject}_freq_transition_main_study.csv').set_index('site_id')
        else:
            freq_transition = pd.read_csv(f'{save_dir}/{subject}_freq_transition_full_study.csv').set_index('site_id')

        info_subject = info.loc[freq_transition.index.values, :]
        freq_transition_info = pd.concat([freq_transition, info_subject],axis = 1)
        freq_transition_info_genes =  freq_transition_info.loc[freq_transition_info['locus_type']  =='CDS',:]
        num_genes_hit[i] = len(freq_transition_info_genes)
        num_snps_change[i] = len(freq_transition_info)
        
        if main_study:
            freq_transition_info.to_csv(f'{save_dir}/{subject}_info_freq_transition_main_study.csv')
        else:
            freq_transition_info.to_csv(f'{save_dir}/{subject}_info_freq_transition_full_study.csv')


        info_subject = info.loc[freq_transition.index.values, :]
        freq_transition_info = pd.concat([freq_transition, info_subject],axis = 1)
        freq_transition_info_genes =  freq_transition_info.loc[freq_transition_info['locus_type']  =='CDS',:]
        num_genes_hit[i] = len(freq_transition_info_genes)
        num_snps_change[i] = len(freq_transition)
        if main_study:
            freq_transition_info.to_csv(f'{save_dir}/{subject}_info_freq_transition_main_study.csv')
        else:
            freq_transition_info.to_csv(f'{save_dir}/{subject}_info_freq_transition_full_study.csv')



    transition_df = pd.DataFrame(data = {'Subject': subjects, 'Num Genes Hit': num_genes_hit,
                                    'Changing sites': num_snps_change})
    transition_df['Species'] = species
    if main_study:
        transition_df.to_csv(f'{save_dir}/gene_info_main_study.csv')
    else:
        transition_df.to_csv(f'{save_dir}/gene_info_full_study.csv')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('out_dir', action='store',
                    help='Outdir prefix where files are stored')
    parser.add_argument('save_dir', action = 'store',
                       help = 'location where to save site_info')
    parser.add_argument('species', action = 'store',
                       help = 'species to perform analysis on')
    parser.add_argument('--main_study',
                    action='store_true')

    args = parser.parse_args()
    species_dir = f'{args.out_dir}/{args.species}'
    save_dir = f'{args.save_dir}/{args.species}'
    if  args.main_study:
        subject_fnames = glob(f'{save_dir}/*freq_transition_main_study.csv')
    else:
        subject_fnames = glob(f'{save_dir}/*freq_transition_full_study.csv')

    subjects = [ fname.split('/')[-1].split('_')[0] for fname in subject_fnames]

    get_main(subjects, species_dir, save_dir, args.species,  main_study = args.main_study)


        
    


