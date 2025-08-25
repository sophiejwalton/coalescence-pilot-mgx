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

def export_plot_pdf(p, fname_prefix):
    p.output_backend = "svg"
    plot_fname_svg=f'plots/{fname_prefix}.svg'
    plot_fname_pdf=f'plots/{fname_prefix}.pdf'
    bokeh.io.export_svg(p,filename=plot_fname_svg)
    convert_scripts = f'rsvg-convert -f pdf -o {plot_fname_pdf} {plot_fname_svg}' 
    os.system(convert_scripts)


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
    
  #  tidy_data['inoculumn'] = tidy_data["sample"].transform(lambda x: e003_metadata.loc[x, 'inoculumn'])
   # tidy_data['inoculumn_sample'] = tidy_data['inoculumn'].transform(get_in)
    
    
    return tidy_data

def make_mesocosm_timecourse(tidy_data, title = '',
                            color = 'grey', alpha = 1.,
                                              limits = (-.5,7.5)   ):
    hv_curve = hv.Curve(data = tidy_data.sort_values(by='passage'),
                kdims=['passage', ],
                vdims=['freq','site_id']
                ).groupby('site_id'
                ).opts(width=500, height=250,
                color = color,
                ylabel = 'Strain AA Allele Frequency',
                       xlabel='Passage',
                title = title,
                #show_grid=True,
                line_width = 1.,
                       alpha = alpha,
               xlim = limits,
                ylim = (-0.05, 1.03)).overlay()
    
    return hv_curve

def get_bootstrap_sel_coeffs(metadata,freq_children, depth_med, depth_children, reads, readsopp, parent_snps1,parent_snps2, inoculumn,species, thresh = 1e-3, n_bootstraps = 1000):
    #med = freq_children.loc[parent_snps,:].median(axis = 0)
    #freq_masked = freq_children.mask((freq_children > freq_thresh * med),axis = 0)
   # freq_masked = freq_masked .mask((freq_masked  < med / freq_thresh),axis = 0)

    mesocosms = []
    initial_freqs=[]
    metadata = metadata.loc[metadata['sample'].isin(freq_children.columns.values),:]
    metadata_non_zero = metadata.loc[metadata['passage']>0,:].sort_values(by='passage')

    for mesocosm in metadata['mesocosm'].unique():
        samples = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'sample'].values)
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values)
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        if inoculumn not in metadata['sample'].values:
            continue 
        random_inds_all = np.random.choice(freq_children.index.values, 1000)
        snps_all_plot = freq_children.loc[random_inds_all, samples]
        tidy=get_tidy_df(snps_all_plot, metadata, value_name = 'freq')
        p1 = make_mesocosm_timecourse(tidy,color=bokeh.palettes.Set2[8][2],alpha = .1)

        random_inds_p1 = parent_snps1.copy()
        if len(parent_snps1) > 1e3:
            random_inds_p1 = np.random.choice(parent_snps1, 1000)

        snps_all_plot = freq_children.loc[random_inds_p1, samples]
        tidy=get_tidy_df(snps_all_plot, metadata, value_name = 'freq')
        p2 = make_mesocosm_timecourse(tidy,color=bokeh.palettes.Set2[8][3],alpha = .1)

        random_inds_p2 = parent_snps2.copy()
        if len(parent_snps2) > 1e3:
            random_inds_p2 = np.random.choice(parent_snps2, 1000)
        
        snps_all_plot = freq_children.loc[random_inds_p2, samples]
        tidy=get_tidy_df(snps_all_plot, metadata, value_name = 'freq')
        p3 = make_mesocosm_timecourse(tidy,color=bokeh.palettes.Set2[8][3],alpha = .1)
        p = hv.render(p1*p2*p3)
        export_plot_pdf(p, f'workflow/report/plots/{species}_{inoculumn}')
    


def get_main(metadata,species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered):
  #  info, depth, freq = load_and_sort_files(species_dir, species)
   # med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
   # good_samples = med_nonzero_depth[med_nonzero_depth>10.]
   # depth = depth[good_samples.index.values]
   # freq = freq[good_samples.index.values]  
   # depth_filtered= depth_filtering(depth)
   # freq_filtered = freq_masked(freq, depth_filtered)

    freq_inoculumns = freq_filtered[parent_samples]
    print(parent_samples)
    
    if len(freq_inoculumns.columns.values)<2:
        return pd.DataFrame(), pd.DataFrame()
    # get distinguishing SNPs for inoculumns - there should be like 1k distinguishing SNPs 
    # use only Alt Allele as marker... so sites where strain allele is alt allele in one strain and not other strain
    # is the marker 
    distinguishing_snps1, distinguishing_snps2= get_distinguishing_snps(freq_inoculumns, thresh = .8)
   # print(len(distinguishing_snps))
    #distinguishing_snps.to_csv('distinguishing_snps.csv')
    parent1_snps = distinguishing_snps1[distinguishing_snps1 == 2].index.values
    parent2_snps = distinguishing_snps2[distinguishing_snps2 == 2].index.values
    freq_children = freq_filtered[child_samples]
    depth_children = depth_filtered[child_samples]
    reads_children = freq_children*depth_children
    reads_children_opp = (1-freq_children)*depth_children
    print(reads_children)
   
    print('before', len(parent1_snps), len(parent2_snps))
    parent1_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent1_snps, thresh = .5, sample_thresh=.75)
    parent2_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent2_snps, thresh = .5, sample_thresh=.75)
    print('after', len(parent1_snps), len(parent2_snps))
    med_depth_children = depth_children.copy().replace(0, np.nan).median(skipna=True)

    count_parent1 = pd.DataFrame(get_count(freq_children, parent1_snps)).rename(columns = {0: parent_samples[0]})  
    count_parent2 = pd.DataFrame(get_count(freq_children, parent2_snps)).rename(columns = {0: parent_samples[1]})  

    ps = get_bootstrap_sel_coeffs(metadata,freq_children, med_depth_children, depth_children, reads_children, reads_children_opp, parent1_snps,parent2_snps, inoculumn,species, n_bootstraps = 1000)
  #  parent2_info = get_bootstrap_sel_coeffs(metadata,freq_children, med_depth_children, reads_children, reads_children_opp, parent2_snps, n_bootstraps = 1000)

   # return pd.concat([count_parent1, count_parent2],axis=1), parent1_info, #parent2_info 


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


    if not path.isdir(save_dir):
        mkdir(save_dir) 
    info, depth, freq = load_and_sort_files(species_dir, args.species)
    #print(info.columns.values)
    #print(info.index.values)
    freq = repolarize_against_reference(freq, info)
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4_good.csv')

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
        'AA-AE-mBHI', 'AA-AF-mBHI',
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    for inoculumn in inoculumn_list:
        print(inoculumn)
        parent_samples, child_samples = get_parent_children(inoculumn,metadata)
      #  print(parent_samples)
#        print(freq_filtered.columns.values)
#AAAA
        parent_samples = list(np.intersect1d(parent_samples, freq_filtered_in.columns.values))
        child_samples = list(np.intersect1d(child_samples, freq_filtered_in.columns.values))
        #child_samples = freq_filtered_in.columns.values
        print(parent_samples)
        print(child_samples)
        if len(parent_samples) < 2:
            continue
        if parent_samples[1] in child_samples:
            child_samples.remove(parent_samples[1])
        if parent_samples[0] in child_samples:
            child_samples.remove(parent_samples[0])
        if parent_samples[0] == 'A2-e003Coalescence-Inoculumn-mBHI':
            parent_samples = ['A2-e003Coalescence-mBHI-inoculumn-redo', parent_samples[1]]
        
        get_main(metadata,species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples])
       
        
    

