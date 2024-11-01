import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
#from tqdm import tqdm
import argparse

from snp_analysis_tools_sherlock import *
import warnings
warnings.filterwarnings('ignore')

def get_species_abundances(species, good_samples):
    species_abundances = pd.read_csv('~/git/household-transmission-mgx/workflow/out/midasOutput/species/species_profile_all.txt', sep = '\t')

    species_abundances_species = species_abundances.loc[species_abundances['species_id'] == species, :]
    species_abundances_samples = species_abundances_species.loc[species_abundances_species['sample'].isin(good_samples), :]
    species_abundances_samples['timepoint'] = species_abundances_samples['sample'].str.split("-").str[-1]
    species_abundances_samples['timepoint'] = pd.to_numeric(species_abundances_samples['timepoint'])
    species_abundances_samples = species_abundances_samples.sort_values('timepoint')
 
    return species_abundances_samples 



  
def get_times_and_ordered_values(median_depth_series,species_abundances_samples,diversity_samples, good_samples):
   # print(median_depth_series)
    median_depth_series = median_depth_series.sort_values('timepoint')

    median_depths = median_depth_series[0].values
    times = median_depth_series['timepoint'].values
    species_abundances_samples = species_abundances_samples.sort_values('timepoint')

    rel_abundances = species_abundances_samples['relative_abundance'].values
    diversity_samples = diversity_samples.sort_values('timepoint')
    diversities = diversity_samples[0].values

    return times, median_depths, rel_abundances, diversities 



def get_good_samples(samples, subject, freq, depth_filtered, badsamples, depth_threshsold = 20, repolarize = True):
    subject_samples = []
    for sample in samples:
        if subject in sample:
            subject_samples.append(sample)
    
    depth_filtered_subject = depth_filtered[subject_samples]
    depth_filtered_subject_nonzero = depth_filtered_subject.copy().replace(0, np.nan)
    depth_filtered_subject_nonzero_median = depth_filtered_subject_nonzero.median()
    depth_filtered_subject_nonzero_median_good = depth_filtered_subject_nonzero_median[depth_filtered_subject_nonzero_median >= depth_threshsold]
   # print('depth filtered')
    good_samples = depth_filtered_subject_nonzero_median_good.index.values
    
    good_timepoints = []
    for sample in good_samples:
        if sample in badsamples:
            good_samples.remove(sample)
        else:
            good_timepoints.append(int(sample.split('-')[-1]))

    depth_filtered_subject_nonzero_median_good = depth_filtered_subject_nonzero_median_good[good_samples]
    freq_polarized = freq.copy()
    if repolarize:
        first_good_time = np.min(np.array(good_timepoints))
        sample_to_polarize = f'HouseholdTransmission-Stool-{subject}-{str(first_good_time).zfill(3)}'
       # print('Polarize')
        freq_polarized = polarize_species(freq_polarized, sample_to_polarize)
    good_freq = freq_polarized[good_samples]
    good_depth = depth_filtered[good_samples]
    return good_freq, good_depth, depth_filtered_subject_nonzero_median_good, good_samples

        
def get_subject_dfs(freq, depth, subject, samples):
    subject_samples = []
    for sample in samples:
        if subject in sample:
            subject_samples.append(sample)
    subject_depth = depth[subject_samples]
    subject_freq = freq[subject_samples]
    
    return subject_freq, subject_depth
def get_diversity_series(freq, depth_filtered, genome_length, sites_considered):
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

def cleanup_and_polarize(freq, median_depth_series, diversity_series, subject):
    good_samples = list(freq.columns.values)
    good_timepoints = []
    ok_timepoints = []
    times = []

    for sample in good_samples:
        if median_depth_series[sample]>= 20:
            if diversity_series[sample] < 1e-3:
                good_timepoints.append(int(sample.split('-')[-1]))
            else:
                ok_timepoints.append(int(sample.split('-')[-1]))
        else: 
          #  print('cool')
            times.append(int(sample.split('-')[-1]))
  #  print(len(good_samples), len(good_timepoints), len(ok_timepoints), len(times))
    if len(good_timepoints) > 0:
        first_good_time = np.min(np.array(good_timepoints))
    elif len(ok_timepoints)> 0:
        first_good_time = np.min(np.array(ok_timepoints))
    else: 
        first_good_time = np.min(times)
        
   # print('YAY', good_samples)
    sample_to_polarize = f'HouseholdTransmission-Stool-{subject}-{str(first_good_time).zfill(3)}'
   # print(sample_to_polarize)
    freq_polarized = polarize_species(freq[good_samples].copy(), sample_to_polarize)
    return freq_polarized, good_samples


def filter_sites_across_samples(good_depth, good_freq):
    
    good_samples = good_depth.columns
    counts = good_depth.count(axis = 1)
    
    mintimes = round(len(good_samples)*.8)
    passing_sites = counts[counts > mintimes]
       # print(len(good_freq))
    good_freq = good_freq.loc[passing_sites.index.values, :]
    good_depth = good_depth.loc[passing_sites.index.values, :]
       # print(len(good_freq))
    return good_depth, good_freq

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


def filter_transition_frequency(moving_snps,  depth_transition, median_depth_series):
    good_sites = []
    for site in moving_snps.index.values:
        good_changing_site = False 
      
        
        coverage_ratio = depth_transition.loc[site, :]/median_depth_series[depth_transition.columns]
        coverage_ratio = coverage_ratio
        site_freq = moving_snps.loc[site, :]
        low_samples  = site_freq[site_freq < .2]
        high_samples = site_freq[site_freq > .8]
        high_samples_coverage = coverage_ratio[site_freq > .8]
        low_samples_coverage = coverage_ratio[site_freq < .2]
        if high_samples_coverage.max()/low_samples_coverage.min() < 2:
            good_changing_site = True 

        else:

            for sample in low_samples.index:
                good_highs = high_samples_coverage[high_samples_coverage < coverage_ratio[sample]*2]

                if len(good_highs) > 0:
                    good_changing_site = True

                    break
        if good_changing_site:
                good_sites.append(site)
    
    freq_filtered_transition = moving_snps.loc[good_sites, :]
    
    return  freq_filtered_transition 

def subsample_and_plot(good_freq, good_depth, color = 'grey',  x_limits = (-10, 70)):
    good_freq_subsampled = good_freq.copy()
    good_depth_subsampled = good_depth.copy()
    if len(good_freq) > 1000:
        good_freq_subsampled, good_depth_subsampled= get_plotting_snps(good_freq, good_depth)
    tidy_data_freq = get_tidy_df(good_freq_subsampled)
    tidy_data_depth = get_tidy_df(good_depth_subsampled, value_name = 'depth')
    tidy_data_freq_good = tidy_data_freq.loc[~np.isnan(tidy_data_depth['depth']), :]
    snps_plot = make_subject_timecourse(tidy_data_freq_good.sort_values('timepoint'),
                                                    color = color,
                                                     limits = x_limits)
    return snps_plot
    
def get_main(all_ids, species_dir, save_dir, subjects, badsamples, species, genome_length):
    info, depth, freq = load_and_sort_files(species_dir)
    info = info.set_index('site_id')
    freq = freq.set_index('site_id')
    depth = depth.set_index('site_id')
    sites_considered = int(len(depth))
    depth_filtered = depth_filtering(depth)
    samples = freq.columns.values
    good_sites = np.zeros(len(subjects))
    plots = []
    abundance_plots = []
    div_plots = []
    depth_plots = []
    snps_plots = []
    titles = []
  
    ids = []

    changing_sites = np.zeros(len(subjects))
    good_changing_sites = np.zeros(len(subjects))
 
    print(subjects)
    for i, subject in enumerate(subjects):
        if ('Y' in subject) or ('Z' in subject):
            continue 

  
       # print("HSIT")
        
       # print(subject)
        abx = False
        if subject[-1] == 'A':
            abx = True
       # if subject != 'XHC':
        #    continue

        subject_freq, subject_depth = get_subject_dfs(freq, depth_filtered, subject, samples)

        diversity_series, _= get_diversity_series(subject_freq, subject_depth, genome_length, sites_considered)
        median_depth_series = subject_depth.copy().replace(0, np.nan).median(skipna = True)
        good_depth_samples = median_depth_series[median_depth_series > 10].index.values
        good_depth_samples = list(good_depth_samples)
    
 
        good_depth_samples_good = []
        for sample in good_depth_samples:
            if sample not in badsamples:
          
                good_depth_samples_good.append(sample)
                
     
      #  print('v2', good_depth_samples_good)

        if len(good_depth_samples_good) == 0:
            continue 
        good_freq, good_samples = cleanup_and_polarize(subject_freq[good_depth_samples_good], median_depth_series, diversity_series, subject)
        
        good_depth = subject_depth[good_samples]
        good_freq = good_freq[good_samples]
        good_depth_na = good_depth*good_depth.isna() + 1. 
        good_freq = good_freq*good_depth_na
        
        good_depth, good_freq = filter_sites_across_samples(good_depth, good_freq)
        
        diversity_series, _ = get_diversity_series(good_freq, good_depth, genome_length, sites_considered)
        good_sites[i] = len(good_depth)
        ### Now grab the transitioning SNPs
        freq_polarized_transition = get_transition_frequency_snps(good_freq, good_depth)
        times = np.array([int(sample.split('-')[-1]) for sample in freq_polarized_transition.columns.values])
        x_range = (-10,np.max(times))
  
        if len(good_freq) == 0:
            continue 
        snps_plot = subsample_and_plot(good_freq, good_depth, color = 'grey', x_limits = x_range)

        #except:
            #continue
        if len(freq_polarized_transition)> 0:
            depth_transition = good_depth.loc[freq_polarized_transition.index.values, :]
 
            freq_transition_filter = filter_transition_frequency(freq_polarized_transition.copy(),  depth_transition, median_depth_series)
            if len(freq_transition_filter) > 0:

                moving_snps = freq_transition_filter.index.values 
                good_stuff = info.loc[moving_snps].drop(columns = ['count_samples', 'count_a', 'count_c', 'count_g', 'count_t'])
                good_stuff = good_stuff.sort_values('site_id')
          #  if len(good_stuff) < 1000:
           #     good_stuff.to_csv(f'{save_dir}/{species}_{subject}_moving_snps.csv')
                depth_transition = depth_transition.loc[freq_transition_filter.index.values, :]
           # depth_transition.to_csv(f'{save_dir}/{species}_{subject}_depth_moving_snps.csv')
            #freq_transition_filter.to_csv(f'{save_dir}/{species}_{subject}_freq_moving_snps.csv')
            
                snps_plot2 = subsample_and_plot(freq_transition_filter, depth_transition, color = 'blue',
                                            x_limits = x_range)
            
                snps_plot = snps_plot*snps_plot2
            
                changing_sites[i] = len(depth_transition)
                good_changing_sites[i] = len(moving_snps)
            
            # now do the whole thing with the later timepoints 
                first_good_time = 64
                times = np.array([int(sample.split('-')[-1]) for sample in freq_transition_filter.columns.values])
                good_times = times[times > first_good_time]
      
                if len(good_times) > 1:

                    followup_samples = [f'HouseholdTransmission-Stool-{subject}-{str(time).zfill(3)}' for time in good_times]

                    good_times_transition = freq_transition_filter[followup_samples]
                    depth_good_times = depth_transition[followup_samples]
                    depth_good_times = depth_good_times.loc[good_times_transition.index.values, :]
    # find new transition 
                    new_transition = get_transition_frequency_snps(good_times_transition, depth_good_times)
                    filtered_good_times = filter_transition_frequency(new_transition,  
                                                      depth_good_times.loc[new_transition.index.values, :],
                                                      median_depth_series[followup_samples])
                    moving_snps = filtered_good_times.index.values
                    good_stuff = info.loc[moving_snps].drop(columns = ['count_samples', 'count_a', 'count_c', 'count_g', 'count_t'])
                    good_stuff = good_stuff.sort_values('site_id')

                  #  good_stuff.to_csv(f'{save_dir}/{species}_{subject}_follow_up_sampling_moving_snps.csv')
                   # depth_follow.to_csv(f'{save_dir}/{species}_{subject}_follow_up_sampling_depth_moving_snps.csv')
                   # filtered_good_times.to_csv(f'{save_dir}/{species}_{subject}_follow_up_sampling_freq_moving_snps.csv')

                    snps_plot2 = subsample_and_plot(freq_transition_filter.loc[moving_snps, :], 

                                                    depth_transition.loc[moving_snps, :], color = 'red', x_limits = x_range)
                    snps_plot = snps_plot*snps_plot2


            
            
            
            
        else:
            changing_sites[i] = 0
            good_changing_sites[i] = 0
            

        median_depth_series = median_depth_series[good_freq.columns.values]
        median_depth_series = pd.DataFrame(median_depth_series).reset_index()
        
        
        median_depth_series['timepoint'] = median_depth_series['index'].str.split("-").str[-1]
        median_depth_series['timepoint'] = pd.to_numeric(median_depth_series['timepoint'])
        species_abundances_samples  = get_species_abundances(species, good_samples)
        diversity_series = pd.DataFrame(diversity_series).reset_index()
        diversity_series['timepoint'] = diversity_series['index'].str.split("-").str[-1]
        diversity_series['timepoint'] = pd.to_numeric(diversity_series['timepoint'])
            
        times, median_depths, rel_abundances, diversities = get_times_and_ordered_values(median_depth_series,
                                                                               species_abundances_samples, diversity_series,
                
                                                                                         good_samples)
       # print(times, median_depths, rel_abundances, diversities)
        #print(len(times), len(median_depths), len(rel_abundances), len(diversities))
        
        depth_plot, abundance_plot, div_plot = depth_abund_div(info, times, median_depths, 
                                                               rel_abundances,diversities, 
                                                               limits= x_range, abx = abx)
        
      #  depth_plot, abundance_plot = depth_abund(info, times, median_depths, 
           #                                          rel_abundances, limits= x_range, abx = abx)

        snps_plot = hv.render(snps_plot)
        #if abx:
        snps_plot.add_layout(bokeh.models.BoxAnnotation(bottom=-0.2, top=0, left = 22, right = 29,
                                                    fill_alpha=1, fill_color='palegreen', level='underlay'))
        snps_plot.add_layout(bokeh.models.BoxAnnotation(bottom=-0.2, top=0, left = 29, right = 34,
                                                    fill_alpha=1, fill_color='#FFD700', level='underlay'))
        snps_plot.add_layout(bokeh.models.BoxAnnotation(bottom=-0.2, top=0, left = 34, right = 40,
                                                    fill_alpha=1, fill_color='lightpink', level='underlay'))
        snps_plot.legend.visible = False
       # snps_plot.title.text = f'{species} in {subject}'
     
        str_ints = ['1', '2', '3', '4']
        titles.append(f'{species} in {subject}')
        new_id = ''.join(list(np.random.choice(str_ints, 4)))
        while new_id in all_ids:
            new_id = ''.join(list(np.random.choice(str_ints, 4)))
        ids.append(new_id)
        all_ids.append(new_id)
        snps_plot.title.text = new_id
        depth_plot.title.text = new_id
        abundance_plot.title.text = new_id
        div_plot.title.text = new_id
           # bokeh.io.export_png(bokeh.layouts.gridplot([depth_plot, abundance_plot, div_plot, snps_plot], ncols = 1),
            #                                           filename = f'{save_dir}/{species}_{subject}.png')
 
            
            
            
        plots.append(bokeh.layouts.gridplot([depth_plot, abundance_plot, div_plot, snps_plot], ncols = 1))
     
       

   # print(plots)
    fraction_filtered_sites = good_sites/genome_length
    subject_info = pd.DataFrame(data = {'Subject': subjects, 'Good Sites': good_sites, 'Changing Sites': changing_sites, 'Good Changing Sites': good_changing_sites, 
                                       'Fraction Filtlered Site': fraction_filtered_sites})
   # subject_info.to_csv(f'{save_dir}/{species}_subject_snps_info.csv')



    titles_ids_df = pd.DataFrame(data = {'ids': ids, 'titles': titles})
    return plots,  titles_ids_df, all_ids
    
  
  #  bokeh.io.export_png(bokeh.layouts.gridplot([depth_plots, abundance_plots, div_plots, snps_plots]), 
   #                     filename = f'{save_dir}/{species}_find_evo_follow.png')
    






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('out_dir', action='store',
                    help='Outdir prefix where files are stored')
    parser.add_argument('save_dir', action = 'store', 
                       help = 'location where to save site_info')

    
    args = parser.parse_args()
    species_info = pd.read_csv('~/git/household-transmission-mgx/workflow/analysis/species_info_full.csv')

    subject_species_info_full = pd.read_csv('~/git/household-transmission-mgx/workflow/analysis/qp_samples_per_species_subject.csv')
    

    

    sampleBlacklist = pd.read_csv('~/git/household-transmission-mgx/workflow/analysis/sampleBlacklist.txt',
                                  sep = ' ')
    contaminated_samples = pd.read_csv('~/git/household-transmission-mgx/workflow/analysis/contaminatedSamples.txt',
                                  sep = ' ')
    badsamples = np.concatenate([sampleBlacklist['sample'].values, contaminated_samples['sample'].values])
    badsamples = [f'HouseholdTransmission-Stool-{sample}' for sample in badsamples]
   # print(badsamples)
    all_plots = []

    all_ids = []
    all_ids = list(pd.read_csv('plots/scrambled_find_evo_follow_panel.csv')['ids'].values)
    all_titles_ids = []
    
    fnames = glob('/Users/lnmerk/git/household-transmission-mgx/workflow/out/midasOutput/snps/HouseholdTransmission-Stool/*/snps_freq.txt.gz')
    og_species_list = ['Bacteroides_ovatus_58035', 'Bacteroides_uniformis_57318', 
                    'Bacteroides_vulgatus_57955', 'Eubacterium_eligens_61678', 'Eubacterium_rectale_56927',
                   'Prevotella_copri_61740']
    
    good_species_info = pd.read_csv('species_subject_info_v2.csv')
    species_list = good_species_info.loc[good_species_info['Subject'] >= 10, 'Species'].values
   # print(len(species_list))
    print(len(species_list))
    #species_list = ['Bacteroides_ovatus_58035', 'Bacteroides_uniformis_57318']
    for i, species in enumerate(species_list):
        print(i, species,)
       # print(species in og_species_list)
       # if species in og_species_list:
         
           # print('DONE')
            #continue
        #species = fname.split('/')[-2]
        
 
            
        subject_species_info  = subject_species_info_full.loc[subject_species_info_full['Species'] == species, :]
        subject_species_info = subject_species_info.loc[subject_species_info['Sample'] > 2, :]
   
        subjects =  subject_species_info['Subject'].values
      
        genome_length = species_info.loc[species_info['species_id'] == species, 'length'].values[0]
        species_dir = f'{args.out_dir}/{species}'
        save_dir = f'{args.save_dir}/{species}'
        plots, titles_ids_df, all_ids = get_main(all_ids, species_dir, save_dir, subjects, badsamples, species, genome_length)
        all_plots = all_plots + plots
    
        all_titles_ids.append(titles_ids_df)
    all_titles_ids_df = pd.concat(all_titles_ids)
    
        
        
    # scramble those plots 
    all_titles_ids_df.to_csv(f'{args.save_dir}/scrambled_find_evo_follow_panel_v2.csv')
       
    np.random.shuffle(all_plots)   
    for i in range(round(len(all_plots)/12)):
        bokeh.io.export_png(bokeh.layouts.gridplot(all_plots[i*12:(i+1)*12], ncols = 4), 
                        filename = f'{args.save_dir}/scrambled_find_evo_follow_panel_{str(i)}_{dec6}.png')
    
    
    
    
    # plot those plots 
    
    
        
    

