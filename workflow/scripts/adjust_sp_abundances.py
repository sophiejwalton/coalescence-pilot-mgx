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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')
    species_df = pd.read_csv('workflow/out/midas2_output/species/species_relative_abundance.tsv',delimiter='\t').set_index('species_id')
    for fname in glob('workflow/out/midas2_output/*/species/species_profile.tsv'):
        sample = fname.split('/')[-3]
        sample_df = pd.read_csv(fname,delimiter='\t').set_index('species_id')
        if sample in species_df.columns.values:
          species_df.loc[sample_df.index.values, sample] = sample_df['marker_relative_abundance'].values	
        #print(sample_df)
       # species_df.loc[sample_df.index.values, sample] = sample_df['marker_relative_abundance'].values 
    species_df.to_csv('workflow/out/midas2_output/species/species_relative_abundance_alt.tsv',sep='\t')

