import pandas as pd
import numpy as np
from os import path

from collections import Counter


from glob import glob
#from tqdm import tqdm



def load_and_sort_files(species_dir,species):
    """
    Takes in snp directory, returns dataframes of various MIDAS2 outputs: info, depth, freq
    """
    info = pd.read_csv(species_dir + f"/{species}.snps_info.tsv", sep = "\t", low_memory=False).set_index('site_id')
    depth = pd.read_csv(species_dir + f"/{species}.snps_depth.tsv", sep = "\t", ).set_index('site_id')
    freq = pd.read_csv(species_dir + f"/{species}.snps_freqs.tsv", sep = "\t", ).set_index('site_id')

    return info, depth, freq

def depth_filtering(depth):
    depth_no_site = depth.copy()
    med = depth_no_site.replace(0, np.nan).median(axis = 0,skipna=True)
    depth_masked_1 = depth_no_site.mask((depth_no_site > 2.5 * med),axis = 0)
    depth_masked = depth_masked_1.mask((depth_masked_1 < med / 2.5),axis = 0)
    depth_masked_absolute = depth_masked.mask(depth_masked < 5)
    return depth_masked_absolute

def freq_masked(freq, depth_filtered):
    depth_filtered_na = depth_filtered*depth_filtered.isna() + 1. 
    return freq.copy()*depth_filtered_na 


def get_diversity_series(freq_filtered, thresh=.2):
    temp_freq = freq_filtered[freq_filtered.notna()].replace(np.nan, .5)
    less_than_thresh = temp_freq < thresh
    less_than_thresh = less_than_thresh.sum()
    greater_than_thresh = temp_freq > 1-thresh
    greater_than_thresh = greater_than_thresh.sum()
    all_nonint_sites = greater_than_thresh + less_than_thresh
    # number of intermediate frequency sites = number total sites with non nan depth - non intermediate frequency sites 
    num_int_sites = freq_filtered.count() - all_nonint_sites
    diversity = num_int_sites/freq_filtered.count() 
    return num_int_sites, diversity


def polarize_species(freq, sample):
    """
    Polarizes all time points to a single sample. If the frequency in the time point to polarize is >0.5,
    all other timepoints at that site are converted to 1-freq.
    Inputted frequency should not have site_id as a column, but as the index.
    """

    samples = list(freq.columns)

    if 'site_id' not in freq.columns:
        freq.index.name = 'site_id'
        # Get site id as a column
        freq.reset_index(inplace = True)

    # Turn into a melted df
    freq_melted = pd.melt(freq, id_vars=['site_id'], value_name = 'freq', var_name = 'sample')
#    print(freq_melted.head())
    freq_melted_polarized = freq_melted.copy()
    sites_to_flip = freq.loc[freq[sample] > 0.5]['site_id']

    freq_melted_polarized['freq'].where(~freq_melted['site_id'].isin(sites_to_flip), 1 - freq_melted['freq'], inplace = True)

    # Return to rectangle dataframe
    freq_polarized = freq_melted_polarized.set_index(['site_id', 'sample'])['freq'].unstack()
  #  print(freq_polarized.head())
    if 'site_id'  in samples:
        samples.remove('site_id')

    # Reorder columns
  #  print(freq_polarized.columns)
#    freq_polarized = freq_polarized[samples]
    
    return freq_polarized


def get_qp_sites(freq, depth_filtered, genome_length, sites_considered):
    temp_freq = freq[depth_filtered.notna()].replace(np.nan, .5)

    less_than_20 = temp_freq < .2
    less_than_20 = less_than_20.sum()

    greater_than_80 = temp_freq > .8
    greater_than_80 = greater_than_80.sum()

    all_nonint_sites = greater_than_80 + less_than_20 
    # number of intermediate frequency sites = number total sites with non nan depth - non intermediate frequency sites 
    num_int_sites = depth_filtered.count() - all_nonint_sites

    # number of bad sites = number sites considered - number good sites 
    bad_sites = sites_considered - depth_filtered.count()

    # total sites = genome length - bad sites 
    total_sites = -bad_sites + genome_length
    qp_sites = num_int_sites/total_sites
    return qp_sites, num_int_sites 
  

def get_important_statistics(depth_filtered, freq, genome_length):
    qp_sites, num_int_sites = get_diversity(freq, depth_filtered, genome_length)

    
    depth_filtered_nonzero = depth_filtered.copy().replace(0, np.nan)
    depth_median = depth_filtered_nonzero.median()
    depth_median = depth.rename('Median Nonzero Depth')
    qp_sites  = qp_sites.rename('Diversity')
    num_int_sites = num_int_sites.rename('Polymorphic Sites')
    
    site_info = pd.concat([qp_sites, num_int_sites, depth_median])
    site_info = pd.DataFrame(data = [qp_sites, num_int_sites, depth]).transpose().reset_index()
    site_info = site_info.rename(columns = {'index': 'Sample'})
    site_info['Subject'] = site_info['Sample'].str.split('-').str[-2]
    site_info['Timepoint'] = site_info['Sample'].str.split('-').str[-1]
    site_info['Household'] = site_info['Subject'].str[1]
    site_info['Study Arm'] = site_info['Subject'].str[0]
    site_info['ABX'] = site_info['Subject'].str[2] == 'A'
    return site_info 
    



def get_transition_frequency_snps(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, .5) # so do not help o

    freq_pass_2 = freq_polarized[(temp_freq<0.2).any(axis=1)]
    depth_pass_2 = depth_filtered[(temp_freq<0.2).any(axis=1)]
  #  print(len(freq_pass_2))
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
        
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, .5)

    freq_polarized_transition= freq_pass_2[(temp_freq > 0.8).any(axis=1)]
    #print(freq_polarized_transition)
    return freq_polarized_transition


def repolarize_against_reference(freq, info):
    '''
    Repolarize so that all allele frequencies are the frequency of the alternative allele 
    '''
    info2 = info.reset_index().copy()
    info2['site_id_adjust'] = info2['site_id'] 
    info2 = info2.set_index('site_id_adjust')
    repolarize = info2['ref_allele'] == info2['minor_allele']

    repolarized = repolarize.copy()
    repolarized[:] = False
    
    repolarized_index = repolarize.loc[repolarize == True].index 
   
    repolarized[repolarized_index] = True
    
    freq_polarized = freq.copy().set_index('site_id')
    
    freq_polarized.loc[repolarized,: ] = 1 - freq_polarized.loc[repolarized,: ]
    return freq_polarized,repolarized, repolarize

    

def get_intermediate_frequency_snps(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, -1)

    freq_pass_1 = freq_polarized[~(temp_freq<0.0).any(axis=1)]
    depth_pass_1 = depth_filtered[~(temp_freq<0.0).any(axis=1)]
   
    freq_pass_2 = freq_pass_1[~(freq_pass_1<0.2).all(axis=1)]
    depth_pass_2 = depth_pass_1[~(freq_pass_1<0.2).all(axis=1)]
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
  
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, 1.1)
    freq_polarized_plotting = freq_pass_2[~(temp_freq > 0.8).all(axis=1)]
    return freq_polarized_plotting


def get_plotting_snps(freq_polarized, depth_filtered):
    freq_polarized_plotting =  get_intermediate_frequency_snps(freq_polarized, depth_filtered)
    freq_polarized_subsampled = freq_polarized_plotting.copy()
    depth_filtered_subsampled = depth_filtered.copy()
    if len(freq_polarized_subsampled) > 500:
      
        freq_polarized_subsampled = freq_polarized_plotting.sample(n=500, random_state=10).sort_index()
        depth_filtered_subsampled = depth_filtered.loc[freq_polarized_subsampled.index, :]
    return freq_polarized_subsampled, depth_filtered_subsampled
        
def get_tidy_df(filtered_freq, value_name = 'freq'):
    filtered_freq.index.name = 'site_id'
    tidy_data = filtered_freq.reset_index().melt(id_vars = ['site_id'], var_name = 'sample', value_name = value_name)
    tidy_data['timepoint'] = tidy_data["sample"].str.split("-").str[3]
    tidy_data['Subject'] = tidy_data["sample"].str.split("-").str[2]
    tidy_data['timepoint'] = pd.to_numeric(tidy_data['timepoint'])
    
    return tidy_data

def make_subject_timecourse(tidy_data, abx = True, title = '',
                            color = 'grey',
                                             limits = []):
   # if title == '':
    #    title = tidy_data['Subject'].values[0]
    if len(limits) == 0:
        end = np.max(tidy_data['timepoint'].values)
        limits = (-10, end)
        
    hv_curve = hv.Curve(data = tidy_data,
                kdims=['timepoint', 'freq'],
                vdims=['site_id']
                ).groupby('site_id'
                ).opts(height = 350,
                width = 800,
                color = color,
                ylabel = 'Allele Frequency',
                title = title,
                #show_grid=True,
                line_width = 0.2,
               xlim = limits,
                ylim = (-0.05, 1.03)).overlay()
    
    return hv_curve

def depth_abund(info, times, median_depths, rel_abundances, abx = True, limits=(-10, 700)):

    if limits[1] == (700):
        end = np.max(times) + 10
        limits = (-10, end)

    p1 = bokeh.plotting.figure(height = 400,width = 800, y_axis_label = 'Median Depth',y_axis_type='log',
                               x_range = limits, )
    p1.line(times, median_depths, line_width = 4, color = 'green')
    p1.circle(times, median_depths, size = 10, color = 'green')

    p2 = bokeh.plotting.figure(height = 400,width = 800,y_range=(10**(-4), 2), y_axis_type='log', y_axis_label = 'Relative Abundance', x_range = limits, tools = [])
    p2.line(times, rel_abundances, line_width = 4, color = 'purple')
    p2.circle(times, rel_abundances, size = 10, color = 'purple')

    p1.toolbar_location = None
    p2.toolbar_location = None
    if abx:

        p1.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 36, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
        p2.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 36, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))

    return p1, p2


def depth_abund_divv2(info, times, median_depths, rel_abundances,diversity, abx = True, limits=(-10, 700)):

    if limits[1] == (700):
        end = np.max(times) + 10
        limits = (-10, end)

    p1 = bokeh.plotting.figure(height = 150,width = 800, y_axis_label = 'Median Depth',y_axis_type='log',
                               x_range = limits, )
    p1.line(times, median_depths, line_width = 4, color = 'green')
    p1.circle(times, median_depths, size = 10, color = 'green')

    p2 = bokeh.plotting.figure(height = 150,width = 800, y_range=(10**(-4), 2), y_axis_type='log', y_axis_label = 'Relative Abundance', x_range = limits, tools = [])
    p2.line(times, rel_abundances, line_width = 4, color = 'purple')
    p2.circle(times, rel_abundances, size = 10, color = 'purple')
    
    
    p3 = bokeh.plotting.figure(height = 150,width = 800, y_range=(10**(-6), 10**(-1)), y_axis_type='log', y_axis_label = 'Diversity', x_range = limits, tools = [])
    p3.line(times, diversity, line_width = 4, color = 'teal')
    p3.circle(times, diversity, size = 10, color = 'teal')

    p1.toolbar_location = None
    p2.toolbar_location = None
    p3.toolbar_location = None
    #if abx:
    

    p1.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p1.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p1.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    
    p2.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p2.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p2.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    
    p3.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p3.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p3.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    return p1, p2, p3


def depth_abund_div(info, times, median_depths, rel_abundances,diversity, abx = True, limits=(-10, 700)):

    if limits[1] == (700):
        end = np.max(times) + 10
        limits = (-10, end)

    p1 = bokeh.plotting.figure(height = 250,width = 800, y_axis_label = 'Median Depth',y_axis_type='log',
                               x_range = limits, )
    p1.line(times, median_depths, line_width = 4, color = 'green')
    p1.circle(times, median_depths, size = 10, color = 'green')

    p2 = bokeh.plotting.figure(height = 250,width = 800, y_range=(10**(-4), 2), y_axis_type='log', y_axis_label = 'Relative Abundance', x_range = limits, tools = [])
    p2.line(times, rel_abundances, line_width = 4, color = 'purple')
    p2.circle(times, rel_abundances, size = 10, color = 'purple')
    
    
    p3 = bokeh.plotting.figure(height = 250,width = 800, y_range=(10**(-6), 10**(-1)), y_axis_type='log', y_axis_label = 'Diversity', x_range = limits, tools = [])
    p3.line(times, diversity, line_width = 4, color = 'teal')
    p3.circle(times, diversity, size = 10, color = 'teal')

    p1.toolbar_location = None
    p2.toolbar_location = None
    p3.toolbar_location = None
    #if abx:
    

    p1.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p1.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p1.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    
    p2.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p2.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p2.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    
    p3.add_layout(bokeh.models.BoxAnnotation(left = 22, right = 29, fill_alpha=1, 
                                                fill_color='palegreen', level='underlay'))
    p3.add_layout(bokeh.models.BoxAnnotation(left = 29, right = 34, fill_alpha=1, 
                                                fill_color='#FFD700', level='underlay'))
    p3.add_layout(bokeh.models.BoxAnnotation(left = 34, right = 40, fill_alpha=1, 
                                                fill_color='lightpink', level='underlay'))
    return p1, p2, p3

def make_subject_timecourse_with_transitions(tidy_data, tidy_data_transition,  abx = True, title = '',
                                             limits = []
                                             ):
   # if title == '':
    #    title = tidy_data['Subject'].values[0]
    if len(limits) == []:
        end = np.max(tidy_data['timepoint'].values)
        x_limits = (-10, end)
        
    hv_curve = hv.Curve(data = tidy_data,
                kdims=['timepoint', 'freq'],
                vdims=['site_id']
                ).groupby('site_id'
                ).opts(height = 350,
                width = 800,
                color = 'grey',
                ylabel = 'Allele Frequency',
                title = title,
                #show_grid=True,
                line_width = 0.2,
               xlim = limits,
                ylim = (-0.05, 1.03)).overlay()
    hv_curve2 = hv.Curve(data = tidy_data_transition,
                kdims=['timepoint', 'freq'],
                vdims=['site_id']
                ).groupby('site_id'
                ).opts(height = 350,
                width = 800,
                color = 'blue',
                ylabel = 'Allele Frequency',
                title = title,
                #show_grid=True,
                line_width = 0.3,
                xlim = limits,
                ylim = (-0.05, 1.03)).overlay()
    p = hv.render(hv_curve*hv_curve2)
    if abx:
        p.add_layout(bokeh.models.BoxAnnotation(bottom=-0.2, top=0, left = 29, right = 36,
                                                fill_alpha=1, fill_color='#FFD700', level='underlay'))
    p.legend.visible = False
    return p

def get_intermediate_frequency_snps_v2(freq_polarized, depth_filtered):
    # only get intermediate frequency snps for plotting
    # Replace nan depth sites with -1, so when you check for all < 0.2, they don't help or detract
    try:
        freq_polarized = freq_polarized.set_index('site_id').copy()
    except:
        
        freq_polarized = freq_polarized.copy()
    
    temp_freq = freq_polarized[depth_filtered.notna()].replace(np.nan, -1)

    freq_pass_2 = freq_polarized[~(temp_freq<0.2).all(axis=1)]
    depth_pass_2 = depth_filtered[~(temp_freq<0.2).all(axis=1)]
                # Replace nan depth sites with 1.1, so when you check for all > 0.8, they don't help or detract
  
    temp_freq = freq_pass_2[depth_pass_2.notna()].replace(np.nan, 1.1)

    freq_polarized_plotting = freq_pass_2[~(temp_freq > 0.8).all(axis=1)]
    return freq_polarized_plotting
