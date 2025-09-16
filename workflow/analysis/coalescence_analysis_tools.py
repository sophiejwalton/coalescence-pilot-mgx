import pandas as pd
import numpy as np

from snp_analysis_tools_sherlock import *

import iqplot
import bokeh.plotting
import bokeh.io
import holoviews as hv
from holoviews import dim, opts
import bokeh.models
from bokeh.layouts import gridplot

hv.extension('bokeh')

def transform_df(df_abundance):
    df_abundance['Lineage'] = df_abundance['species_id'].transform(lambda x: df_metadata.loc[df_metadata['species_id'] == x,'Lineage'].values[0])
    df_abundance['species'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-1])
    df_abundance['genus'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-2])
    df_abundance['family'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-3])
    df_abundance['phyla'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[1])
    return df_abundance


def get_inoculumn_sort(x):
    subjects = list(np.sort(x.split('-')[:-1]))
    media =  x.split('-')[-1]
    return '-'.join(subjects + [media])



    
def analyze_fitness(diversity_df1,minor_strain, major_strain,
                               minor_strain_subject, major_strain_subject):
    new_df = []
    for i, type_meso in enumerate(diversity_df1['type_meso'].unique()):
        mesos = diversity_df1.loc[diversity_df1['type_meso'] ==type_meso, 'mesocosm'].unique()
        df_type_meso = diversity_df1.loc[diversity_df1['type_meso'] == type_meso,:]
        for mesocosm in mesos:
            
            df_meso = diversity_df1.loc[diversity_df1['mesocosm'] == mesocosm,:]
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


def analyze_diversity(diversity_df1,minor_strain, minor_strain_subject, separate_plots_per_mesocosm=True,
                     real_species = ''):
    
    e003_metadata = pd.read_csv('e003_metadata_cultures_round2.csv').drop(columns = 'Unnamed: 0')
    e003_metadata['type_meso'] = e003_metadata['type_mesocosm']


    in_df = e003_metadata.loc[e003_metadata['is_inoculumn'],:].set_index('inoculumn').copy()

    in_series = in_df['sample']
    in_dict = in_series.to_dict()
    p2 = bokeh.plotting.figure(width = 500, height = 300)
    palette = bokeh.palettes.Set2[8]
    for i, type_meso in enumerate(diversity_df1['type_meso'].unique()):
        mesos = diversity_df1.loc[diversity_df1['type_meso'] ==type_meso, 'mesocosm'].unique()
        
     #   p3 = bokeh.plotting.figure(width = 300, height = 200)
        df_type_meso = diversity_df1.loc[diversity_df1['type_meso'] == type_meso,:]
        color = palette[0]
        if type_meso.split('-')[-1] == 'mBHI':
            color = palette[1]
        p2.circle(df_type_meso['passage'].values,
                   df_type_meso[minor_strain].values,color = color, legend_label = type_meso, size = 5)
        for mesocosm in mesos:
            
            df_meso = diversity_df1.loc[diversity_df1['mesocosm'] == mesocosm,:]
            in_sample = df_meso['inoculumn_sample'].unique()[0]
            #print(type_meso, in_sample)
            if in_sample in diversity_df1.index.values:
                df_meso.loc[in_sample,:] = diversity_df1.loc[in_sample,:]
           # df_meso = pd.concat([df_meso, diversity_df1.loc[diversity_df1['sample'] == in_sample,:]])
            df_meso = df_meso.sort_values(by = 'passage')
          #  p.line(df_meso['passage'].values,
           #        df_meso[minor_strain].values, legend_label)

                #print(np.max(df_meso['passage'].values))
               # print(len(df_meso))
                
            p2.line(df_meso['passage'].values,
                   df_meso[minor_strain].values,color=color)
                

              
    p2.yaxis.axis_label = 'Minor strain abundance'
    p2.xaxis.axis_label = 'Passage'
    inoculumn_stuff = '-'.join(type_meso.split('-')[:-1])
    p2.title.text = f'{real_species} {minor_strain_subject} in {inoculumn_stuff}'
    p2.y_range = bokeh.models.Range1d(-.005,1.005)
    p2.x_range = bokeh.models.Range1d(-.5,7.5)

      #  plots.append(p2)

            

  #  if separate_plots_per_mesocosm:
    return p2
        
    