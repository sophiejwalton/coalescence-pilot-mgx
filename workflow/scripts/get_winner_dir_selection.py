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
    
def transform_df(df_abundance):
    df_abundance['Lineage'] = df_abundance['species_id'].transform(lambda x: df_metadata.loc[df_metadata['species_id'] == x,'Lineage'].values[0])
    df_abundance['species'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-1])
    df_abundance['genus'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-2])
    df_abundance['family'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-3])
    df_abundance['phyla'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[1])
    return df_abundance


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
    distinguishing_snps = get_distinguishing_snps(freq_inoculumns, thresh = .8)
    #print(distinguishing_snps)
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


def analyze_fitness(df_info,minor_strain, major_strain,
                               minor_strain_subject, major_strain_subject):
    new_df = []
    for i, type_meso in enumerate(df_info['type_meso'].unique()):
        mesos = df_info.loc[df_info['type_meso'] ==type_meso, 'mesocosm'].unique()
        df_type_meso = df_info.loc[df_info['type_meso'] == type_meso,:]
        for mesocosm in mesos:
            
            df_meso = df_info.loc[df_info['mesocosm'] == mesocosm,:]
            in_sample = df_meso['inoculumn_sample'].unique()[0]
            df_meso['shift_from_inoculumn'] = np.nan
           # df_meso[f'shift {minor_strain_subject}'] = np.nan
            df_meso[f'opp_strain_shift_from_inoculumn'] = np.nan
            if in_sample in df_info.index.values:
                df_meso.loc[in_sample,:] = df_info.loc[in_sample,:]
                df_meso['shift_from_inoculumn'] = df_meso[minor_strain] - df_info.loc[in_sample,minor_strain]
            #    df_meso[f'shift {minor_strain_subject}'] = df_meso[minor_strain] - df_info.loc[in_sample,minor_strain]
                df_meso[f'opp_strain_shift_from_inoculumn'] = df_meso[major_strain] - df_info.loc[in_sample,major_strain]
                
        
            df_meso = df_meso.sort_values(by = 'passage')
            new_df.append(df_meso)
            
    
    return pd.concat(new_df)

def get_inoculumn_sort(x):
    subjects = list(np.sort(x.split('-')[:-1]))
    media =  x.split('-')[-1]
    return '-'.join(subjects + [media])

def get_in(x):
    #print('yay',x)
    if x in list(in_dict.keys()):
     
        return in_dict[x]
    else:
        return ''

def get_species_tracking(species):
    fnames = glob(f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snps/{species}/*_parent_freqs.csv')
    
    all_dfs = []
    for info_fname in np.sort(fnames): 
       # print(info_fname)
       # fname = f'~/git/coalescence-pilot-mgx/workflow/report/track_snps/{species}/AA-AE-mGAM_parent_freqs.csv'
        df = pd.read_csv(info_fname).rename(columns = {'Unnamed: 0': 'sample'})
        # print(df)
        if len(df) ==0:
            continue

        medium_df = df.set_index('sample').median()
        minor_strain = medium_df.index.values[0]
        major_strain = medium_df.index.values[1]
        # strains = df.columns.values[1:]
        #  subjects = [e003_metadata.loc[e003_metadata['sample'] == strain,'parent_subjects'].values[0].split('-')[0] for strain in strains]
        # print(subjects)

        df_info = pd.concat([df.set_index('sample'), 
                                e003_metadata.loc[e003_metadata['sample'].isin(df['sample'].unique()),:].set_index('sample')],axis=1)
        df_info['inoculumn_sample'] = df_info['inoculumn'].transform(get_in)
            #print(minor_strain)
        minor_strain_subject =  e003_metadata.loc[e003_metadata['sample'] == minor_strain,
                                              'parent_subjects'].values[0].split('-')[0]
        major_strain_subject =  e003_metadata.loc[e003_metadata['sample'] == major_strain,
                                              'parent_subjects'].values[0].split('-')[0]


        # print(minor_strain)
    #    p2 = analyze_diversity(df_info, minor_strain, 
     #                           e003_metadata.loc[e003_metadata['sample'] == minor_strain,'parent_subjects'].values[0].split('-')[0])
        df2 = analyze_fitness(df_info, minor_strain, major_strain,
                               minor_strain_subject, major_strain_subject)
        all_dfs.append(df2[['experiment',
        'sub_experiment', 'comm', 'parent_subjects', 'is_inoculumn',
        'parent_media', 'media', 'passage', 'mesocosm', 'single-subject',
        'inoculumn', 'type_meso', 'inoculumn_sample',  'opp_strain_shift_from_inoculumn',
        'shift_from_inoculumn']].reset_index())
       # bokeh.io.show(p2)
      #  plots_all.append(p2)

    if len(all_dfs) == 0:
        return 
    df_full = pd.concat(all_dfs)
  #  df_full_p5p7 = df_full.loc[df_full['passage'] > 3,:].sort_values(by = 'sample')
    df_full['dir_of_shift'] = 1*(df_full['shift_from_inoculumn'] > 0.) + -1*(df_full['shift_from_inoculumn'] < 0.)
    df_full['zero_change'] = df_full['shift_from_inoculumn'] == 0.
    df_full['total_shift'] = np.abs(df_full['shift_from_inoculumn'] + df_full['opp_strain_shift_from_inoculumn'])
    df_full.to_csv(f'{save_dir}/{species}_shifts.csv', index = False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
   # parser.add_argument('--species', action = 'store', 
    #                   help = 'species to perform analysis on')
#    parser.add_argument('--inoculumn', action = 'store', 
 #                      help = 'inoculumn')
    args = parser.parse_args()
   # species_dir = f'{args.indir}/{args.species}'
   # save_dir = f'{args.outdir}/{args.species}'
    save_dir = args.outdir

#    parent_samples, child_samples = get_parent_children(args.inoculumn)


    if not path.isdir(save_dir):
        mkdir(save_dir) 

    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    fnames = glob(f'{args.indir}/*/')

    species_list = [fname.split('/')[-2] for fname in fnames]
    e003_metadata = pd.read_csv('/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/analysis/e003_coal_metadata_full.csv').drop(columns = 'Unnamed: 0')

    e003_metadata['inoculumn'] = e003_metadata['inoculumn'].transform(get_inoculumn_sort)

    in_df = e003_metadata.loc[e003_metadata['is_inoculumn'],:].set_index('inoculumn').copy()

    in_series = in_df['sample']
    in_dict = in_series.to_dict()
 
    for species in species_list:
        print(species)
        get_species_tracking(species)





       
        
    

