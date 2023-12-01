import glob 
import pandas as pd
import argparse



if __name__ == '__main__':
    fname = '/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midas2_output/merge/species/species_relative_abundance.tsv'
    df_rel_abundance = pd.read_csv(fname, delimiter = '\t')
    df_rel_abundance_melted = pd.melt(df_rel_abundance, var_name = 'Sample', id_vars='species_id', value_name = 'relative_abundance')
    abundant_species_df = df_rel_abundance_melted.groupby(['species_id']).max().reset_index()
    abundant_species_df =  abundant_species_df.sort_values(by = 'relative_abundance', ascending = False)
    abundant_species_df = abundant_species_df[:150]
    abundant_species_df.to_csv('/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midas2_output/abundant_species.csv')

   