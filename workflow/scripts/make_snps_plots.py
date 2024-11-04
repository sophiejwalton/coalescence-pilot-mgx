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


def get_distinguishing_snps(freq_inoculumns, thresh = .8):
    # find snps that are greater than .8 (aka alternative alleles > .8)
    detect_df = freq_inoculumns > thresh
    detect_df['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df = detect_df.loc[detect_df['site_present'] == 0,:]
   # print(detect_df[freq_inoculumns.columns.values[0]])
    detect_df['diff'] = detect_df[freq_inoculumns.columns.values[0]].astype(int) - detect_df[freq_inoculumns.columns.values[1]].astype(int)
    return detect_df['diff']

def get_frequency_parent(freq_children, parent_snps):
    median_freq = freq_children.loc[parent_snps,:].median(axis = 0)

    return median_freq 


def get_parent_children(inoculumn):
#    metadata = pd.read_csv('config/e003_coal_metadata_full.csv')
    metadata = pd.read_csv('config/e003_metadata_cultures_round2.csv')
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
        return pd.DataFrame(), [], []
    # get distinguishing SNPs for inoculumns - there should be like 1k distinguishing SNPs 
    # use only Alt Allele as marker... so sites where strain allele is alt allele in one strain and not other strain
    # is the marker 
    distinguishing_snps = get_distinguishing_snps(freq_inoculumns, thresh = .999)
    print(distinguishing_snps)
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
  #  freq_parent1.to_csv('freq_parent1.csv')
  #  freq_parent1['parent'] = parent_samples[0]
    freq_parent2 = get_frequency_parent(freq_children, parent2_snps)
    freq_parent2 = pd.DataFrame(freq_parent2).rename(columns = {0: parent_samples[1]})



# freq_parent2 = freq_parent2.T
   # freq_parent2['parent'] = parent_samples[1]
   # print(freq_parent2)
    return pd.concat([freq_parent1, freq_parent2],axis=1), parent1_snps, parent2_snps


def subsample_and_plot(good_freq, good_depth, color = 'grey',  x_limits = (-.5, 7.5), alpha = 1.):
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
    
    tidy_data['inoculumn'] = tidy_data["sample"].transform(lambda x: e003_metadata.loc[x, 'inoculumn'])
    tidy_data['inoculumn_sample'] = tidy_data['inoculumn'].transform(get_in)
    
    
    return tidy_data
def get_inoculumn_sort(x):
    subjects = list(np.sort(x.split('-')[:-1]))
    media =  x.split('-')[-1]
    return '-'.join(subjects + [media])





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
    e003_metadata = pd.read_csv('workflow/analysis/e003_coal_metadata_full.csv').drop(columns = 'Unnamed: 0')

    e003_metadata['inoculumn'] = e003_metadata['inoculumn'].transform(get_inoculumn_sort)

    in_df = e003_metadata.loc[e003_metadata['is_inoculumn'],:].set_index('inoculumn').copy()

    in_series = in_df['sample']
    in_dict = in_series.to_dict()
    e003_metadata = e003_metadata.set_index('sample')
    def get_in(x):
    #print('yay',x)
    	if x in list(in_dict.keys()):
     
       		return in_dict[x]
    	else:
        	return ''
    if not path.isdir(save_dir):
        mkdir(save_dir) 
    info, depth, freq = load_and_sort_files(species_dir, args.species)
    freq = repolarize_against_reference(freq, info)
    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>10.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    plots = []
    for inoculumn in inoculumn_list:
        print(inoculumn)
        parent_samples, child_samples = get_parent_children(inoculumn)
        parent_samples = list(np.intersect1d(parent_samples, freq_filtered.columns.values))
        child_samples = list(np.intersect1d(child_samples, freq_filtered.columns.values))
        _, freq_filtered_in = filter_sites_across_samples(depth_filtered[parent_samples + child_samples], freq_filtered[parent_samples + child_samples].copy(), )

        freq_parents, parent1_snps, parent2_snps = get_main(species_dir,args.species, parent_samples, child_samples, freq_filtered_in)
       # if len(freq_parents) == 0:
          
        inoculumnstr = ''.join(inoculumn.split('/'))
        freq_parents.to_csv(f'{save_dir}/{inoculumnstr}_parent_freqs.csv')
        print(len(parent1_snps), 'party')
        
        mesocosms = e003_metadata.loc[e003_metadata['inoculumn'] == inoculumn, 'mesocosm'].unique()
#        print(mesocosms, 'MESOCOSM')
        plots = []
        if len(parent1_snps) >1000:
            parent1_snps = np.random.choice(parent1_snps, 1000)
        if len(parent2_snps)>1000:
            parent2_snps = np.random.choice(parent2_snps, 1000)
        for i, mesocosm in enumerate(mesocosms):
            print(mesocosm)
           # if len(freq_parents) == 0:
            #   continue 
          #  mesocosm = ''.join(mesocosm.split('/'))
            samples = e003_metadata.loc[e003_metadata['mesocosm'] == mesocosm, :].index.values
            inoculumn_sample = get_in(inoculumn)
#            print('yay', inoculumn_sample)
 #           print('wee', samples)
            samples = list(samples) + [inoculumn_sample]   
           # print(freq_filtered.columns.values) 
            samples = list(np.intersect1d(samples, freq_filtered.columns.values))     
            freq_filtered_mesocosm = freq_filtered[samples]
            random_snps = np.random.choice(freq_filtered_mesocosm.index.values, 10000)
            freq_filtered_mesocosm_rand  = freq_filtered_mesocosm.loc[random_snps, :]
            freq_filtered_mesocosm_rand = get_tidy_df(freq_filtered_mesocosm_rand, e003_metadata, )
           # print('go', freq_filtered_mesocosm_rand) 
            if i == 1:
                p = make_mesocosm_timecourse(freq_filtered_mesocosm_rand, title = mesocosm )
#            print(p)
                freq_filtered_mesocosm_marker_parent1  = get_tidy_df(freq_filtered_mesocosm.loc[parent1_snps, :], e003_metadata)
                freq_filtered_mesocosm_marker_parent2  = get_tidy_df(freq_filtered_mesocosm.loc[parent2_snps, :], e003_metadata)
                p1 = make_mesocosm_timecourse(freq_filtered_mesocosm_marker_parent1 ,title = f'{mesocosm} parent1',  color = bokeh.palettes.Accent[3][1], alpha = .2 )
                p2 = make_mesocosm_timecourse(freq_filtered_mesocosm_marker_parent2 , color = bokeh.palettes.Accent[3][2], title = f'{mesocosm} parent2', alpha = .2 )
                mesocosmstr = ''.join(mesocosm.split('/'))
                bokeh.io.export_png(bokeh.layouts.gridplot([hv.render(p), hv.render(p1), hv.render(p2)],ncols = 1), filename = f'{save_dir}/{inoculumnstr}_{mesocosm}_snps.png')
            
           # plots.append([hv.render(p), hv.render(p1), hv.render(p2)])
           # print(plot, 'YAY')
       # if len(plots)>0:
#        bokeh.io.export_png(bokeh.layouts.gridplot(plots,ncols = 1), filename = f'{save_dir}/{inoculumn}_snps.png')










    


       
        
    

