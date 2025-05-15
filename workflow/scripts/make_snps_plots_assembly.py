import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
import itertools as it
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
import warnings
#import iqplot
import bokeh.plotting
import bokeh.io
import holoviews as hv
from holoviews import dim, opts
import bokeh.models
from bokeh.layouts import gridplot

hv.extension('bokeh')
warnings.filterwarnings('ignore')


def get_parent_children_assembly(inoculumn,metadata):
#    metadata = pd.read_csv('config/e003_coal_metadata_full.csv')
 #  metadata = pd.read_csv('config/e003_metadata_cultures_round2.csv')
    child_samples = list(metadata.loc[metadata['inoculumn'] == inoculumn, 'sample'].values)
    

    parent_subjects = inoculumn.split('-')[:-1]
    parent_media = inoculumn.split('-')[-1]
    ins = metadata.loc[metadata['is_inoculumn'],:]
    ins = ins.loc[ins['parent_media'] == parent_media,:]
    print(parent_subjects)
    ins1 = ins.loc[ins['parent_subjects'] == parent_subjects[0] + '-' +  parent_subjects[0],:]

    ins2 = ins.loc[ins['parent_subjects'] == parent_subjects[1] + '-' +  parent_subjects[1],:]
    if (len(ins1) ==0) or (len(ins2) ==0): 
        return np.nan, np.nan
    parent_samples = [ins1['sample'].values[0], ins2['sample'].values[0]]
    media = inoculumn.split('-')[-1]
    in_parent1 = f'{parent_subjects[0]}-{parent_subjects[0]}-{media}'
    in_parent2 = f'{parent_subjects[1]}-{parent_subjects[1]}-{media}'
    child_samples_parent1_ss =  list(metadata.loc[metadata['inoculumn'] == in_parent1, 'sample'].values)
    child_samples_parent2_ss =  list(metadata.loc[metadata['inoculumn'] == in_parent2, 'sample'].values)
    return parent_samples, child_samples # + child_samples_parent1_ss + child_samples_parent2_ss 

def subsample_and_plot(good_freq, good_depth, color = 'grey',  x_limits = (-.5, 5.5), alpha = 1.):
    good_freq_subsampled = good_freq.copy()
    good_depth_subsampled = good_depth.copy()
    if len(good_freq) > 1000:
        good_freq_subsampled, good_depth_subsampled= get_plotting_snps(good_freq, good_depth)
    tidy_data_freq = get_tidy_df(good_freq_subsampled)
    tidy_data_depth = get_tidy_df(good_depth_subsampled, value_name = 'depth')
    tidy_data_freq_good = tidy_data_freq.loc[~np.isnan(tidy_data_depth['depth']), :]
    snps_plot = make_mesocosm_timecourse(tidy_data_freq_good.sort_values('passage'),
                                                    color = color,
                                        alpha = alpha,
                                                     limits = x_limits)
    return snps_plot
    
def get_tidy_df(filtered_freq, e003_metadata, value_name = 'freq'):
    filtered_freq.index.name = 'site_id'
    tidy_data = filtered_freq.reset_index().melt(id_vars = ['site_id'], var_name = 'sample', value_name = value_name)
    tidy_data = tidy_data.loc[tidy_data['sample'].isin(e003_metadata.index.values)]

    tidy_data['passage'] = tidy_data["sample"].transform(lambda x: e003_metadata.loc[x, 'passage'])
    tidy_data['passage'] = pd.to_numeric(tidy_data['passage'])

    tidy_data['mesocosm'] = tidy_data["sample"].transform(lambda x: e003_metadata.loc[x, 'mesocosm'])
    
    
    return tidy_data

def get_moving(freq_filtered_mesocosm ):
    df_counts = freq_filtered_mesocosm.count(axis=1)
    fixed_at_zero = freq_filtered_mesocosm==0.
    fixed_at_zero = fixed_at_zero.sum(axis=1)
    fixed_at_zero = fixed_at_zero[fixed_at_zero < df_counts].index.values


    fixed_at_one=freq_filtered_mesocosm==1.
    fixed_at_one = fixed_at_one.sum(axis=1)
    fixed_at_one= fixed_at_one[fixed_at_one < df_counts].index.values

    good_sites = np.intersect1d(fixed_at_zero, fixed_at_one)

    return freq_filtered_mesocosm.loc[good_sites,:]

def make_mesocosm_timecourse(tidy_data, title = '',
                            color = 'grey', alpha = 1.,
                                              limits = (-.5,7.5)   ):
   # if title == '':
    #    title = tidy_data['Subject'].values[0]
   # if len(limits) == 0:
    #    end = np.max(tidy_data['timepoint'].values)
     #   limits = (-10, end)
    
   
    hv_curve = hv.Curve(data = tidy_data,
                kdims=['passage', 'freq'],
                vdims=['site_id']
                ).groupby('site_id'
                ).opts(height = 350,
                width = 800,
                color = color,
                ylabel = 'Allele Frequency',
                title = title,
                #show_grid=True,
                line_width = 0.2,
                       alpha = alpha,
               xlim = limits,
                ylim = (-0.05, 1.03)).overlay()
    
    return hv_curve

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
    assembly_metadata = pd.read_csv('workflow/analysis/assembly_glycerol_metadata_with_redo.csv').drop(columns='Unnamed: 0')#.set_index('sample')
    assembly_metadata['pseudo_time_passage'] = assembly_metadata['passage'] 
    assembly_metadata.loc[assembly_metadata['pseudo_time_passage'] ==1,'pseudo_time_passage'] = 6. 
    assembly_metadata=assembly_metadata.loc[assembly_metadata['sample']!='Assembly-G12-fecal-AA-fecal-0_S703',:]

    if not path.isdir(save_dir):
        mkdir(save_dir) 
    info, depth, freq = load_and_sort_files(species_dir, args.species)
    freq = repolarize_against_reference(freq, info)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh=2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    assembly_metadata=assembly_metadata.loc[assembly_metadata['pseudo_time_passage']<6,:]

    _, freq_filtered_in = filter_sites_across_samples(depth_filtered, freq_filtered.copy(), thresh=.75)
    plots = []
    for mesocosm in assembly_metadata['mesocosm'].unique():
        subject = mesocosm.split('-')[1]
        fecal_meso= f'fecal-{subject}-fecal'
        samples = assembly_metadata.loc[assembly_metadata['mesocosm'] == mesocosm, 'sample'].values
        sample_fecal= assembly_metadata.loc[assembly_metadata['mesocosm'] == fecal_meso, 'sample'].values[0]
        samples = list(samples) + [sample_fecal]   
           # print(freq_filtered.columns.values) 
        print(samples)
        print(freq_filtered_in)
        samples = list(np.intersect1d(samples, list(freq_filtered_in.columns.values)))     
        freq_filtered_mesocosm = freq_filtered_in[samples]
        freq_filtered_mesocosm =get_moving(freq_filtered_mesocosm )
        if sample_fecal in freq_filtered_mesocosm.columns.values:
            freq_filtered_mesocosm = polarize_species(freq_filtered_mesocosm, sample_fecal)
        freq = repolarize_against_reference(freq, info)
        if len(freq_filtered_mesocosm)==0:
            continue

        random_snps = np.random.choice(freq_filtered_mesocosm.index.values, 10000)
        freq_filtered_mesocosm_rand  = freq_filtered_mesocosm.loc[random_snps, :]
        freq_filtered_mesocosm_rand = get_tidy_df(freq_filtered_mesocosm_rand, assembly_metadata.set_index('sample'))
           # print('go', freq_filtered_mesocosm_rand) 
        if len(mesocosm)>0:
            mesocosm = ''.join(mesocosm.split('/'))
            p = hv.render(make_mesocosm_timecourse(freq_filtered_mesocosm_rand, title = f'{args.species} in {mesocosm}' ))
            bokeh.io.export_png(p, filename = f'{save_dir}/{args.species}_{mesocosm}_assembly_snps.png')
            plots.append(p)
    bokeh.io.export_png(bokeh.layouts.gridplot(plots, ncols=2), filename = f'{save_dir}/{args.species}_assembly_snps.png')









    


       
        
    

