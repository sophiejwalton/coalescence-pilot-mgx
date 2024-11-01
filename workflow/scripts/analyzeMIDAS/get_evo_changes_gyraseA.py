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

def is_gyraseA(x):
    return 'GO:0003918' in x


def get_gyraseA_genes(genome_info):
    genome_info['isGyraseA']  = genome_info['functions'].transform(is_gyraseA)
    gyraseAgenes = genome_info.loc[genome_info['isGyraseA'], 'gene_id'].valuess
    
    return gyraseAgenes


 
def get_main(subjects, species_dir, save_dir, species,  gyraseAgenes=[], main_study = True):
    info, _, _ = load_and_sort_files(species_dir)
    info = info.set_index('site_id')

    good_changing_sites = np.zeros(len(subjects))
    gyraseAgenes = np.zeros((len(subjects), len(gyraseAsites)))
    if len(gyraseAsites) == 0:
        data = {'Subject': [], 'Good Changing Sites': []}
        if main_study:
            data.to_csv(f'{save_dir}/gyraseA_df_main_study.csv')
        else:
            data.to_csv(f'{save_dir}/gyraseA_df_full_study.csv')
        return None


    for i, subject in enumerate(subjects):

        if args.main_study:
            freq_transition_info = pd.read_csv(f'{save_dir}/{subject}_info_freq_transition_main_study.csv').set_index('site_id')
        else:
            freq_transition_info = pd.read_csv(f'{save_dir}/{subject}_info_freq_transition_full_study.csv').set_index('site_id')

        good_changing_sites[i] = len(freq_transition_info)
        for j, gene  in  enumerate(gyraseAgenes):
            if gene  in freq_transition_info['gene_id'].values:
                gyraseAgenes[i][j] = 1

        freq_transition_info_gyrase = freq_transition_info.loc[freq_transition_info['gene_id'].isin(gyraseAgenes), :]
        if len(freq_transition_info_gyrase) > 0:
            if args.main_study:
                freq_transition_info_gyrase.to_csv(f'{save_dir}/{subject}_gyraseA_freq_transition_main_study.csv')
            else:
                freq_transition_info_gyrase.to_csv(f'{save_dir}/{subject}_gyraseA_freq_transition_full_study.csv')

    data = {'Subject': subjects, 'Good Changing Sites': good_changing_sites}
    for i, site in enumerate(gyraseAsites):
        data[str(site)] = gyraseA_sites[:,i].flatten()
    gyraseA_df = pd.DataFrame(data = data)
    gyraseA_df['Species'] = species
    if main_study:
        gyraseA_df.to_csv(f'{save_dir}/gyraseA_df_main_study.csv')
    else:
        gyraseA_df.to_csv(f'{save_dir}/gyraseA_df_full_study.csv')



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

    subject_fnames = glob(f'{save_dir}/*good_freq.csv.gz')
    subjects = [ fname.split('/')[-1].split('_')[0] for fname in subject_fnames]

    genome_info = pd.read_csv('/home/groups/bhgood/midas_db_v1.2/rep_genomes/{args.species}/genome.features.gz', compression  =  'gzip',delimeter ='\t')
    gyraseAgenes = get_gyraseA_genes(genome_info)


    get_main(subjects, species_dir, save_dir, args.species,  gyraseAgenes = gyraseAgenes, main_study = args.main_study)
    


        
    

