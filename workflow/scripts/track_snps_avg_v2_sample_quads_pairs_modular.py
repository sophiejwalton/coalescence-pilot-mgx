import pandas as pd
import numpy as np
from os import path, mkdir
from glob import glob
import argparse
import itertools as it
from scipy.stats import linregress
from snp_analysis_tools_sherlock import *
from evo_changes_tools import *
from track_snps_funcs import *
import warnings
warnings.filterwarnings('ignore')




def get_bootstrap_sel_coeffs(metadata, freq_children, depth_med, parent_snps, n_bootstraps = 10000, initial_sample = 1, final_sample=3):
    #med = freq_children.loc[parent_snps,:].median(axis = 0)
    #freq_masked = freq_children.mask((freq_children > freq_thresh * med),axis = 0)
   # freq_masked = freq_masked .mask((freq_masked  < med / freq_thresh),axis = 0)
    snps = freq_children.loc[parent_snps,:]
    print(snps.columns.values)
    boot_med = []
    boot_low = []
    boot_high = []
    act_med = []
    samples_good = []
    metadata = metadata.loc[metadata['sample'].isin(snps.columns.values),:]
    metadata_non_zero = metadata.loc[metadata['passage']>0,:].sort_values(by='passage')
    n_snps = len(parent_snps)

    diff7_info_lower=[]
    diff7_info_upper=[]
    diff7_info_med = []

    diff5_info_lower=[]
    diff5_info_upper=[]
    diff5_info_med = []

    diffboth_info_lower=[]
    diffboth_info_upper=[]
    diffboth_info_med = []

    pred_7_upper = []
    pred_7_lower= []
    pred_5_upper = []
    pred_5_lower= []

    pred7_med = []
    pred5_med = []
    mesocosms = []
    as_extreme7 = []
    as_extreme5 = []
    act_7 = []
    act_5 = []
    for mesocosm in metadata['mesocosm'].unique():
        meso_df = metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, :]
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values) + [0]
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        #if inoculumn not in metadata['sample'].values:
         #   continue 
        if final_sample not in passages:
            continue 
        if initial_sample==1:
            sample1 = meso_df.loc[meso_df['passage']==1,'sample'].values
            dt=final_sample - initial_sample
            dt5=4
            dt7=6
        elif initial_sample == 0:
            sample1 = inoculumn 
            dt=final_sample - initial_sample
            dt5=5
            dt7=7
        else: # both!!!!! 
            sample0 = inoculumn 
            sample1 = meso_df.loc[meso_df['passage']==1,'sample'].values
            dt=2
            dt5=4
            dt7=6
      #  print(sample1)      
        if sample1 not in metadata['sample'].values:
            continue 
#        if initial_sample == 'both' and (inoculumn not in metadata['sample'].values):
 #           continue 

        sample3 = meso_df.loc[meso_df['passage']==final_sample,'sample'].values
        p5_act = np.nan
        p7_act = np.nan
        if 5 in passages:
            sample5 = meso_df.loc[meso_df['passage']==5,'sample'].values
            p5_act = get_freq_est(snps.loc[parent_snps,sample5].values, depth_med[sample5].values[0])
        if 7 in passages:
            sample7 = meso_df.loc[meso_df['passage']==7,'sample'].values
            p7_act = get_freq_est(snps.loc[parent_snps,sample7].values, depth_med[sample7].values[0])
        
        diffs7 = np.zeros(n_bootstraps)
        diffs5 = np.zeros(n_bootstraps)
        diff_both = np.zeros(n_bootstraps)

        pred_7_ests = np.zeros(n_bootstraps)
        pred_5_ests = np.zeros(n_bootstraps)
        est_7s_ests = np.zeros(n_bootstraps)
        est_5s_ests = np.zeros(n_bootstraps)
        
        for sample in range(n_bootstraps):
            bs_sample = np.random.choice(parent_snps, size=len(parent_snps),replace=False)
            int1 = round(len(bs_sample)/4)
            int2 = round(2*len(bs_sample)/4)
            int3 = round(3*len(bs_sample)/4)
            if initial_sample== 'both':
                int1 = round(len(bs_sample)/5)
                int2 = round(2*len(bs_sample)/5)
                int3 = round(3*len(bs_sample)/5)
                int4=round(4*len(bs_sample)/5)

            passage1_est = get_freq_est(snps.loc[bs_sample[:int1],sample1].values, 
                        depth_med[sample1] )
            passage3_est = get_freq_est(snps.loc[bs_sample[int1+1:int2],sample3].values, 
                       depth_med[sample3].values[0] )
            passage5_est=np.nan
            passage7_est=np.nan
            if 5 in passages:
                passage5_est = get_freq_est(snps.loc[bs_sample[int2+1:int3],sample5].values, 
                        depth_med[sample5].values[0] )
            if 7 in passages:
                passage7_est = get_freq_est(snps.loc[bs_sample[int3+1:],sample7].values, 
                        depth_med[sample7].values[0] )
            
      
            sel_coeff_13 = (1/dt)*(np.log(passage3_est/(1-passage3_est)) + \
                                np.log((1-passage1_est)/(passage1_est)))

            if initial_sample=='both':
                passage1_est = get_freq_est(snps.loc[bs_sample[int4+1:],sample1].values, 
                        depth_med[sample1] )
                passage7_est = get_freq_est(snps.loc[bs_sample[int3+1:int4],sample7].values, 
                        depth_med[sample7].values[0] )

                xs = np.array([0, 1, 3])
                ys = np.array([np.log(passage0_est/(1-passage0_est)),
                             np.log(passage1_est/(1-passage1_est)),
                             np.log(passage3_est/(1-passage3_est))
                             ])

           

                sel_coeff_13 = linregress(xs,ys).slope


            pred_p7 = np.exp(dt7*sel_coeff_13)
            pred_p7 = pred_p7*passage1_est/(1 - passage1_est + passage1_est*pred_p7)

            pred_p5 = np.exp(dt7*sel_coeff_13)
            pred_p5 = pred_p5*passage1_est/(1 - passage1_est + passage1_est*pred_p5)

          #  print(pred_p7, passage7_est)
            diffs7[sample] = pred_p7 - passage7_est 
            diffs5[sample] = pred_p5 - passage5_est
            diff_both[sample] = ((pred_p7 - passage7_est ) + (pred_p5 - passage5_est))/2

            est_7s_ests[sample] = passage7_est
            pred_7_ests[sample] = pred_p7 
            est_5s_ests[sample] = passage5_est
            pred_5_ests[sample] = pred_p5 

        diff7_info_upper.append(np.quantile(diffs7, q = .975))
        diff7_info_lower.append(np.quantile(diffs7, q = .025))

        diff5_info_upper.append(np.quantile(diffs5, q = .975))
        diff5_info_lower.append(np.quantile(diffs5, q = .025))
        diffboth_info_upper.append(np.quantile(diff_both, q = .975))
        diffboth_info_lower.append(np.quantile(diff_both, q = .025))

        pred_7_upper.append(np.quantile(pred_7_ests, q = .975))
        pred_7_lower.append(np.quantile(pred_7_ests, q = .025))
        pred_5_upper.append(np.quantile(pred_5_ests, q = .975))
        pred_5_lower.append(np.quantile(pred_5_ests, q = .025))

        pred7_med.append(np.median(pred_7_ests))
        pred5_med.append(np.median(pred_5_ests))

        diff7_info_med.append(np.median(diffs7))
        diff5_info_med.append(np.median(diffs5))
        diffboth_info_med.append(np.median(diff_both))

        mesocosms.append(mesocosm)
       # p7_act = get_freq_est(snps.loc[parent_snps,sample7].values, depth_med[sample7].values[0])
        act_7.append(p7_act)
        #p5_act = get_freq_est(snps.loc[parent_snps,sample5].values, depth_med[sample5].values[0])
        act_5.append(p5_act)
        as_extreme7.append(np.sum(p7_act>pred_7_ests))
        as_extreme5.append(np.sum(p5_act>pred_5_ests))

    df = pd.DataFrame(data = { 'mesocosms': mesocosms, 'diff7_med': diff7_info_med, 'diff7_lower': diff7_info_lower, 
                                'diff7_upper': diff7_info_upper, 
                                'diff5_med': diff5_info_med, 'diff5_lower': diff5_info_lower, 
                                'diff5_upper': diff5_info_upper, 
                                'pred_5_upper': pred_5_upper, 'act_5': act_5,
                                'pred_5_lower': pred_5_lower, 
                                'pred_7_upper': pred_7_upper, 'act_7': act_7,
                                'pred_7_lower': pred_7_lower, 
                                'pred7_med': pred7_med,'as_extreme7': as_extreme7, 
                                'pred5_med': pred5_med,'as_extreme5': as_extreme5, 
                                'diff_both_lower': diffboth_info_lower, 'diff_both_upper': diffboth_info_upper, 
                                'diff_both_med': diffboth_info_med,
                                })
    df['n_snps'] = n_snps
    df['sample_init'] = initial_sample
    df['final_sample'] = final_sample
    return df
       


def get_main(metadata, species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered,sample_init= 1, final_sample=3):
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
   
    print('before', len(parent1_snps), len(parent2_snps))
    parent1_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent1_snps, thresh = .5, sample_thresh=.75)
    parent2_snps = filter_distinguishing_snps(freq_filtered[child_samples], parent2_snps, thresh = .5, sample_thresh=.75)
    print('after', len(parent1_snps), len(parent2_snps))
    med_depth_children = depth_children.copy().replace(0, np.nan).median(skipna=True)

   # count_parent1 = pd.DataFrame(get_count(freq_children, parent1_snps)).rename(columns = {0: parent_samples[0]})  
   # count_parent2 = pd.DataFrame(get_count(freq_children, parent2_snps)).rename(columns = {0: parent_samples[1]})  

    parent1_info = get_bootstrap_sel_coeffs(metadata, freq_children, med_depth_children, parent1_snps, n_bootstraps = 10000,
                 initial_sample=sample_init, final_sample=final_sample)
    parent2_info = get_bootstrap_sel_coeffs(metadata, freq_children, med_depth_children, parent2_snps, n_bootstraps = 10000,
                initial_sample=sample_init, final_sample=final_sample)
# freq_parent2 = freq_parent2.T
   # freq_parent2['parent'] = parent_samples[1]
   # print(freq_parent2)
    return parent1_info, parent2_info 

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
    metadata = pd.read_csv('workflow/analysis/e003_with_passage_one_redo.csv')
   # metadata = pd.read_csv('workflow/analysis/e003_metadata_cultures_round2_change_AA.csv')
    

    med_nonzero_depth = depth.copy().replace(0, np.nan).median(skipna=True)
    med_nonzero_depth.to_csv(f'{save_dir}/{args.species}_median_depths.csv')
    good_samples = med_nonzero_depth[med_nonzero_depth>=5.]
    depth = depth[good_samples.index.values]
    freq = freq[good_samples.index.values]
    depth_filtered= depth_filtering(depth, depth_thresh = 2.5)
    freq_filtered = freq_masked(freq, depth_filtered)
    inoculumn_list = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
    depth_filtered_in, freq_filtered_in = filter_sites_across_samples(depth_filtered, 
        freq_filtered,thresh=.75)
    sample_init = 1
    final_sample = 3
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
        
        parent1_info, parent2_info = get_main(metadata,species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples],
            sample_init = sample_init, final_sample=final_sample)

        inoculumn = ''.join(inoculumn.split('/'))
    
      #  count_parents.to_csv(f'{save_dir}/{inoculumn}_count_parents.csv')
        parent1_info.to_csv(f'{save_dir}/{inoculumn}_{str(sample_init)}_{str(final_sample)}_parent1_info.csv')
        parent2_info.to_csv(f'{save_dir}/{inoculumn}_{str(sample_init)}_{str(final_sample)}_parent2_info.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

