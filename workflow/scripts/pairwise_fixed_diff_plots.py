import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
import itertools as it
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
from holoviews.operation import histogram
import iqplot
import bokeh.plotting
import bokeh.io
import holoviews as hv
from holoviews import dim, opts
import bokeh.models
from bokeh.layouts import gridplot
from glob import glob

hv.extension('bokeh')
import warnings
warnings.filterwarnings('ignore')


def get_fixed_diffs(freq_filtered,sample):
    freq_filtered_sample = freq_filtered[[sample]].copy()

    good_sites_lower = freq_filtered_sample.loc[freq_filtered_sample[sample] < .2].index.values
    fixed_diffs_lower_snps = freq_filtered.loc[good_sites_lower,:] >.8 
    print(fixed_diffs_lower_snps)
    fixed_diffs_lower_snps = fixed_diffs_lower_snps.sum(axis=0)
    print(fixed_diffs_lower_snps)
    good_sites_upper = freq_filtered_sample.loc[freq_filtered_sample[sample] > .8].index.values
    fixed_diffs_upper_snps = freq_filtered.loc[good_sites_upper,:] <.8 
    fixed_diffs_upper_snps = fixed_diffs_upper_snps.sum(axis=0)

    num_sites_non_int = (freq_filtered >.8).sum(axis=0) + (freq_filtered <.2).sum(axis=0) 
    fixed_diffs = fixed_diffs_upper_snps+fixed_diffs_lower_snps
    fixed_diffs = fixed_diffs.rename('fixed_diffs')
    num_sites_non_int = num_sites_non_int.rename('num_int_sites')
 
    return pd.concat([fixed_diffs,num_sites_non_int],axis=1)


def get_fd_plot(freq_filtered, s1, s2):
    points = hv.Points(freq_filtered, vdims = [s1,s2],kdims=[s1,s2])
    xhist, yhist = (histogram(freq_filtered[s1].values,dimension=dim).opts(logy=True) *
                histogram(freq_filtered[s2].values, dimension=dim).opts(logy=True)
                for dim in [s1,s2])
    composition = (points) << yhist.opts(width=125) << xhist.opts(height=125)

    return points

def get_main(species_dir, save_dir, species,metadata):
    info, depth, freq = load_and_sort_files(species_dir, species)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]

    good_samples = np.intersect1d(good_samples.index.values, metadata['sample'].values)
    depth = depth[good_samples]
    freq = freq[good_samples]  

    depth_filtered= depth_filtering(depth,depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)

    s1s = []
    s2s = []
    all_plots = []
    i=0

    for s1,s2 in it.combinations(freq_filtered.columns.values,2):
        plot = get_fd_plot(freq_filtered[[s1,s2]],s1,s2)

        all_plots.append(hv.render(plot))
        s1s.append(s1)
        s2s.append(s2)
        i=i+1
       # if i >10:
        bokeh.io.export_png(hv.render(plot),
    filename=f'{save_dir}/{species}_{s1}_{s2}_fixed_diffs.png')
        #    break 
  #  ss_df['Strain Shift'] = ss_df['fixed_diffs'] > 1000
  #  if '/' in parent_subjects_media:
   #     parent_subjects_media = ''.join(parent_subjects_media.split('/'))
    bokeh.io.export_png(bokeh.layouts.gridplot(all_plots, ncols = 4), 
    filename=f'{save_dir}/{species}_fixed_diffs.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='basic filtering of sites')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='Outdir prefix where to save stuff')
    parser.add_argument('--indir', action = 'store', 
                       help = 'location where to get stuff from')
    parser.add_argument('--species', action = 'store', 
                       help = 'species to perform analysis on')
    metadata = pd.read_csv('workflow/analysis/e003_with_passage_one_redo.csv')
    metadata = pd.read_csv('workflow/analysis/assembly_glycerol_metadata_with_redo.csv')


    args = parser.parse_args()
    species_dir = f'{args.indir}/{args.species}'
    save_dir = f'{args.outdir}/{args.species}'
    if not path.isdir(save_dir):
        mkdir(save_dir)
    get_main(species_dir, save_dir, args.species, metadata)
    


       
        
    

