import glob 
import pandas as pd
import argparse
from os import path,mkdir
import numpy as np
midasdir = '/oak/stanford/groups/relman/users/kxue/household-transmission-mgx'

def getGoodAbundanceFile(speciesAbundances_small):  
 # speciesAbundances_small = speciesAbundances.loc[speciesAbundances['sample']==sample,:]
  speciesAbundances_small['new_rel_abundance'] = speciesAbundances_small['relative_abundance']/speciesAbundances_small['relative_abundance'].sum()
  speciesAbundances_small['relative_abundance']  = speciesAbundances_small['new_rel_abundance'] 
  return speciesAbundances_small


def get_contig_abundance_file(abundance_df, contig_df):

    species_to_use =list(np.intersect1d(abundance_df['rep_species'].values, contig_df['rep_species'].values)) 

    abundance_df_small = abundance_df.loc[abundance_df['rep_species'].isin(species_to_use),:]
    abundance_df_small['new_rel_abundance'] = abundance_df_small['relative_abundance']/abundance_df_small['relative_abundance'].sum()
    contig_df_small = contig_df.loc[contig_df['rep_species'].isin(species_to_use),:]
   # print('SKX01058' in contig_df['species_id'].unique(), 'woow')
    contig_df_small_gr = contig_df_small.groupby('rep_species').sum()
    contig_df_small['full_genome_length'] = contig_df_small['rep_species'].transform(lambda x: contig_df_small_gr.loc[x,'contig_length'])
    contig_df_small['contig_length_frac'] = contig_df_small['contig_length']/contig_df_small['full_genome_length']
    contig_df_small['genome_rel_abundance'] = contig_df_small['rep_species'].transform(lambda x: abundance_df_small.loc[abundance_df_small['rep_species'] == x, 'new_rel_abundance'].values[0])
    contig_df_small['contig_rel_abundance'] =  contig_df_small['contig_length_frac']*contig_df_small['genome_rel_abundance'] 
    return contig_df_small


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Make genome fasta file')

    # add arguments
    parser.add_argument('--species_info', action='store',
                    help='location of species_list_file')
    #parser.add_argument('--vulgatus_replacement_genome', action='store',)
    parser.add_argument('--outdir', action='store',
                    help='where to store altered abundance file')
    parser.add_argument('--genomedir', action='store',
                    help='where genomes stored at')

    args = parser.parse_args()

    speciesAbundancesog=pd.read_csv(args.species_info)
    speciesAbundances = getGoodAbundanceFile(speciesAbundancesog)
    genome_info = pd.read_csv(f'{args.genomedir}/genomes_info.csv').rename(columns={'species_id':'rep_species'})
    species_info_prev = args.species_info.split('/')[-1].split('.csv')[0]

    theta_genomes = ['GUT_GENOME000472','GUT_GENOME001553','GUT_GENOME034180','GUT_GENOME001329','GUT_GENOME172048',
    'GUT_GENOME001637'] # sort of random selection 
    theta_genomes=['GUT_GENOME034180','GUT_GENOME001329']#'GUT_GENOME000472','GUT_GENOME001553',
                   # 'GUT_GENOME034180','GUT_GENOME001329']
    #theta - GUT_GENOME001120
    # let's start out by only mixing 50:50 with reference genome 
    for theta_genome in theta_genomes:
        for i,genome_abundance in enumerate([0,1e-5,1e-4,1e-3,1e-2,1e-1,.5,.9,1-1e-2,1-1e-3,1-1e-4,1-1e-5, 1]):
            outdir = f'{args.outdir}/{theta_genome}'
            if not path.isdir(outdir):
                mkdir(outdir )
            print(theta_genome, genome_abundance)
            speciesAbundances_small = speciesAbundances.copy()
            theta_abundance = speciesAbundances_small.loc[speciesAbundances_small['rep_species'] == 'GUT_GENOME001120', 'relative_abundance'].values[0]
            speciesAbundances_small.loc[speciesAbundances_small['rep_species'] == 'GUT_GENOME001120', 'relative_abundance'] = theta_abundance*(1-genome_abundance)
            
            other_species = speciesAbundances_small.loc[speciesAbundances_small['rep_species'] != 'GUT_GENOME001120', 'rep_species'].unique()

            df_Bv_small =  speciesAbundances_small.loc[speciesAbundances_small['rep_species'] == 'GUT_GENOME001120', :].copy()
            df_Bv_small['rep_species'] = theta_genome
            df_Bv_small['relative_abundance'] = theta_abundance*genome_abundance
            speciesAbundances_small = speciesAbundances_small._append(df_Bv_small)

            speciesAbundances_small.to_csv(f'{outdir}/{species_info_prev}_{theta_genome}_{str(i)}_abundances.csv')
            contig_info_small= get_contig_abundance_file(speciesAbundances_small, genome_info)
            contig_info_small[['contig_name', 'contig_rel_abundance']].to_csv(f'{outdir}/{species_info_prev}_{theta_genome}_{str(i)}_abundances.txt',sep = '\t',
            header=False,index=False)
          #  contig_info_small.to_csv(f'{outdir}/{args.sample}/Bacteroides_vulgatus_57955_{vulgatus_replacement_genome}_og_contig_abundance.csv',index=False)    
            

