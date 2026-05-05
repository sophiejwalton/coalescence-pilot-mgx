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


def get_parent_snps(freq_inoculumns,thresh=.8):
    # find snps that are greater than .8 (aka alternative alleles > .8)
    detect_df = freq_inoculumns >= thresh
    return detect_df[detect_df].index.values
def get_distinguishing_snps(freq_inoculumns, thresh = .8):
    # find snps that are greater than .8 (aka alternative alleles > .8)
    detect_df = freq_inoculumns >= thresh
    detect_df_opp = freq_inoculumns <= 1 - thresh
    detect_df['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df_opp['site_present'] = freq_inoculumns[freq_inoculumns.columns.values[0]].isna() + freq_inoculumns[freq_inoculumns.columns.values[1]].isna()
    detect_df = detect_df.loc[detect_df['site_present']==0,:]
    detect_df_opp = detect_df_opp.loc[detect_df_opp['site_present']==0,:]
    detect_df['diff1'] = detect_df[freq_inoculumns.columns.values[0]].astype(int) + detect_df_opp[freq_inoculumns.columns.values[1]].astype(int)
    detect_df['diff2'] = detect_df[freq_inoculumns.columns.values[1]].astype(int) + detect_df_opp[freq_inoculumns.columns.values[0]].astype(int)
    return detect_df['diff1'], detect_df['diff2']

def filter_distinguishing_snps(freq_children, parent_snps, thresh = .5, sample_thresh=.75):
    med = freq_children.loc[parent_snps,:].median(axis = 0) 
    diff_from_med = freq_children - med
   # print(med)
    #print('PARTY')
    #print(diff_from_med)
    freq_masked = freq_children.mask((diff_from_med.abs() > thresh),axis=0)

    min_samples = round(len(freq_children.columns.values)*sample_thresh)
    passing_sites = freq_masked.count(axis = 1)
    passing_sites = passing_sites[passing_sites>= min_samples].index.values
    parent_snps_good = np.intersect1d(passing_sites, parent_snps)
    return parent_snps_good

def get_abuns(genome,sample_input,sample_infer, freqs, depths):
    p1 = f'{genome}_{sample_input}_annotated_abundance_{genome}_14_abundances_trial_000'
    p2 = f'{genome}_{sample_input}_annotated_abundance_{genome}_0_abundances_trial_000'
    d1,d2=get_distinguishing_snps(freqs[[p1,p2]])
    d1 = d1[d1==2].index.values
    d2=d2[d2==2].index.values
    d1=filter_distinguishing_snps(freqs, d1)
    estimated_abun_full = []
    estimated_abun_meds = []
    estimated_abun_03 = []
    estimated_abun_34 = []
    reads_vals3 = []
    reads_vals4 = []
    reads_vals0 = []

    inputs = []
    for i,abun in enumerate([0,.0001, .001,.01,.1,.2,.4,.5,.6,.8,.9,.99,.999,.9999, 1., .002,.005,.007,.992,.995,.997,.02,.05,.07,.98,.95,.93,.7]):
        sample=  f'{genome}_{sample_infer}_annotated_abundance_{genome}_{i}_abundances_trial_000'
        vals = freqs.loc[d1,sample]
        reads_vals = ((freqs.loc[d1,sample])*depths.loc[d1,sample])
        reads_opp_vals = ((1-freqs.loc[d1,sample])*depths.loc[d1,sample])
        reads_vals = reads_vals[~np.isnan(vals)].round()
        reads_opp_vals = reads_opp_vals[~np.isnan(vals)].round()
        vals = vals[~np.isnan(vals)]
        est_med= np.median(vals)
        depth_est = depths[sample]
        depth_est_med = depth_est[depth_est>0].median() 
        est = est_med
        est_03 = est_med
        est_34 = est_med
        reads3 = np.nan
        reads4=np.nan
        reads0=np.nan
        if est ==0:
            b = np.sum(reads_vals==3)/np.sum(reads_vals==0)
            est_03 = ((b*6)**(1/3))/depth_est_med
            b = np.sum(reads_vals==4)/np.sum(reads_vals==3)
            est_34 = b*4/depth_est_med
            est=est_03
            reads3 = np.sum(reads_vals==3)
            reads4 = np.sum(reads_vals==4)
            reads0 = np.sum(reads_vals==0)
            if np.sum(reads_vals==3) >100:
                est = est_34
        if est == 1:
            b = np.sum(reads_opp_vals==3)/np.sum(reads_opp_vals==0)
            est_03 = 1-((b*6)**(1/3))/depth_est_med
            est = est_03
            b = np.sum(reads_opp_vals==4)/np.sum(reads_opp_vals==3)
            est_34 = 1-b*4/depth_est_med
            reads3 = np.sum(reads_opp_vals==3)
            reads4 = np.sum(reads_opp_vals==4)
            reads0 = np.sum(reads_opp_vals==0)
            if np.sum(reads_opp_vals==3) >100:
                est = est_34
        print(abun,est)
        estimated_abun_full.append(est)
        estimated_abun_03.append(est_03)
        estimated_abun_34.append(est_34)
        estimated_abun_meds.append(est_med)
        reads_vals3.append(reads3)
        reads_vals4.append(reads4)
        reads_vals0.append(reads0)
        inputs.append(abun)
    df = pd.DataFrame(data={'est':estimated_abun_full, 'input':inputs, 'est_meds': estimated_abun_meds, 'est_03': estimated_abun_03, 'est_34': estimated_abun_34,
                                'reads_vals3': reads_vals3, 'reads_vals4':reads_vals4, 'reads_vals0': reads_vals0})
    df['genome']=genome
    df['sample_input']=sample_input
    df['sample_infer']=sample_infer
    return df #,d1
    

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
    freq = repolarize_against_reference(freq, info)

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)


    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    depth_filtered_in.to_csv(f'{species_dir}/depth_filtered.csv.gz',compression='gzip')
    freq_filtered_in.to_csv(f'{species_dir}/freq_filtered.csv.gz',compression='gzip')
    genomes = ['GUT_GENOME000472','GUT_GENOME001553','GUT_GENOME001637']
    sample_input = 'G8-e003Coalescence-mBHI-p7_S125'
    sample_infer = 'A10-e003Coalescence-mBHI-p7'
    dfs = []
    for sisf in [('G8-e003Coalescence-mBHI-p7_S125','A10-e003Coalescence-mBHI-p7'),
                ('A10-e003Coalescence-mBHI-p7','G8-e003Coalescence-mBHI-p7_S125',)]:
        sample_input,sample_infer =sisf
        for genome in genomes:
            dfsm = get_abuns(genome,sample_input,sample_infer, freq_filtered_in, depth_filtered_in)
            dfs.append(dfsm)
    df = pd.concat(dfs)
    df.to_csv(f'{save_dir}/sims_freqs.csv')


       
        
    

