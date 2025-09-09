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




def get_sample_freq(freqs, reads,readsopp, depth_med):
    reads = reads[~np.isnan(freqs)]
    readsopp = readsopp[~np.isnan(freqs)]
    freqs = freqs[~np.isnan(freqs)]
    freq_est = np.median(freqs)
    if freq_est == 0:
        freq_est = (np.sum(reads==1)/np.sum(reads==0))/depth_med
    if freq_est == 1:
        freq_est = 1-(np.sum(readsopp==1)/np.sum(readsopp==0))/depth_med
    return freq_est 

def get_sample_freq_adjust(freqs,reads,readsopp, sample, snps1,snps2, depth_med):
    est1 = get_sample_freq(freqs.loc[snps1,sample], reads.loc[snps1,sample],
                                    readsopp.loc[snps1,sample], depth_med[sample])
    est2 = get_sample_freq(freqs.loc[snps2,sample], reads.loc[snps2,sample],
                                    readsopp.loc[snps2,sample], depth_med[sample])
    if est1<1e-3:
        est1=1e-3
    if est2<1e-3:
        est2=1e-3
    if est1>1-1e-3:
        est1 = 1-1e-3
    if est2>1-1e-3:
        est2=1-1e-3
    return est1/(est1+est2)



def get_bootstrap_sel_coeffs(metadata,freq_children, depth_med, depth_children, reads_children, reads_children_opp, parent_snps1,parent_snps2, inoculumn, n_bootstraps = 10000, thresh = 1e-3, early_interval = (1,2), late_interval = (6,7)):
    early1,early2 = early_interval
    late1,late2 = late_interval
    metadata = metadata.loc[metadata['sample'].isin(freq_children.columns.values),:]
    metadata_non_zero = metadata.loc[metadata['passage']>0,:].sort_values(by='passage')

    diff7_info_lower=[]
    diff7_info_upper=[]
    diff7_info_med = []

    pred_7_upper = []
    pred_7_lower= []
    pred_7_med = []

    act_7_upper = []
    act_7_lower= []
    act_7_med = []

    early_sel_lower = []
    early_sel_upper = []
    early_sel_med = []

    late_sel_lower = []
    late_sel_upper = []
    late_sel_med = []

    as_extreme7 = [] 
    mesocosms = []

    for mesocosm in metadata['mesocosm'].unique():
        meso_df = metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, :]
        passages = list(metadata_non_zero.loc[metadata_non_zero['mesocosm'] == mesocosm, 'passage'].values) + [0]
        inoculumn = metadata.loc[metadata['mesocosm'] == mesocosm, 'inoculumn_sample'].values[0]
        #if inoculumn not in metadata['sample'].values:
         #   continue 
        if early1 not in passages:
            continue
        if early2 not in passages:
            continue
        if late1 not in passages:
            continue
        if late2 not in passages:
            continue 

        if early1>0:
            sample1 = meso_df.loc[meso_df['passage']==early1,'sample'].values
            dt=early2-early1
            dt7=late2- early2
            dtlate = late2-late1
        elif early2 == 0:
            sample1 = inoculumn 
            dt=early2-early1
            dt7=late2- early2
            dtlate = late2-late1
        sample2 = meso_df.loc[meso_df['passage'] == early2,'sample'].values
        sample3 = meso_df.loc[meso_df['passage'] == late1,'sample'].values
        sample4 = meso_df.loc[meso_df['passage'] == late2,'sample'].values
   
        if sample1 not in metadata['sample'].values:
            continue 

        diffs7 = np.zeros(n_bootstraps)
        pred_7_ests = np.zeros(n_bootstraps)
        est_7s_ests = np.zeros(n_bootstraps)
        late_sels = np.zeros(n_bootstraps)
        early_sels = np.zeros(n_bootstraps)
        
        for sample in range(n_bootstraps):
            bs_sample1 = np.random.choice(parent_snps1, size=len(parent_snps1),replace=False)
            bs_sample2 = np.random.choice(parent_snps2, size=len(parent_snps2),replace=False)

            int11 = round(len(bs_sample1)/4)
            bssample11 = bs_sample1[:int11]

            int21 = round(2*len(bs_sample1)/4)
            bssample21 = bs_sample1[int11+1:int21]
            int31 = round(3*len(bs_sample1)/4)
            bssample31 = bs_sample1[int21+1:int31]
            bssample41 = bs_sample2[int31+1:]

            int12 = round(len(bs_sample2)/4)
            bssample12 = bs_sample1[:int12]
            int22 = round(2*len(bs_sample2)/4)
            bssample22 = bs_sample1[int12+1:int12]
            int32 = round(3*len(bs_sample2)/4)
            bssample32 = bs_sample1[int22+1:int32]
            bssample42 = bs_sample1[int32+1:]

            passage1_est=get_sample_freq_adjust(freq_children,reads_children, reads_children_opp, sample1, bssample11,bssample12, depth_med)

            passage2_est=get_sample_freq_adjust(freq_children,reads_children, reads_children_opp, sample2, bssample21,bssample22, depth_med)

            passage3_est=get_sample_freq_adjust(freq_children,reads_children, reads_children_opp, sample3, bssample31,bssample32, depth_med)
            
            passage4_est=get_sample_freq_adjust(freq_children,reads_children, reads_children_opp, sample4, bssample41,bssample42, depth_med)

            sel_coeff_early = (1/dt)*(np.log(passage2_est/(1-passage2_est)) + \
                                np.log((1-passage1_est)/(passage1_est)))
            sel_coeff_late =  (1/(dtlate))*(np.log(passage4_est/(1-passage4_est)) + \
                                np.log((1-passage3_est)/(passage3_est)))

            pred_p7 = np.exp(dt7*sel_coeff_early)
            pred_p7 = pred_p7*passage2_est/(1 - passage2_est + passage2_est*pred_p7)

          #  print(pred_p7, passage7_est
            
            est_7s_ests[sample] = passage4_est
            pred_7_ests[sample] = pred_p7 
            if pred_p7 <thresh:
                pred_p7=thresh
            if pred_p7 >1-thresh:
                pred_p7=1-thresh
          #  pred_7_ests[sample] = pred_p7 
            diffs7[sample] = pred_p7 - passage4_est 
            late_sels[sample] = sel_coeff_late
            early_sels[sample] = sel_coeff_early

        print(diffs7)

        diff7_info_upper.append(np.quantile(diffs7, q = .975))
        diff7_info_lower.append(np.quantile(diffs7, q = .025))
        diff7_info_med.append(np.median(diffs7))

        pred_7_upper.append(np.quantile(pred_7_ests, q = .975))
        pred_7_lower.append(np.quantile(pred_7_ests, q = .025))
        pred_7_med.append(np.median(pred_7_ests))
        mesocosms.append(mesocosm)
       # p7_act = get_freq_est(snps.loc[parent_snps,sample7].values, depth_med[sample7].values[0])
        act_7_upper.append(np.quantile(est_7s_ests, q = .975))
        act_7_lower.append(np.quantile(est_7s_ests, q = .025))
        act_7_med.append(np.median(est_7s_ests))

        early_sel_upper.append(np.quantile(early_sels, q = .975))
        early_sel_lower.append(np.quantile(early_sels, q = .025))
        early_sel_med.append(np.median(early_sels))

        late_sel_upper.append(np.quantile(late_sels, q = .975))
        late_sel_lower.append(np.quantile(late_sels, q = .025))
        late_sel_med.append(np.median(late_sels))
        #p6_act = get_freq_est(snps.loc[parent_snps,sample6].values, depth_med[sample6].values[0])
        pred_7_ests[pred_7_ests<1e-3]=1e-3
        pred_7_ests[pred_7_ests>1-1e-3]=1-1e-3
        as_extreme7.append(np.sum(est_7s_ests>pred_7_ests))

    df = pd.DataFrame(data = { 'mesocosms': mesocosms, 'diff7_med': diff7_info_med, 'diff7_lower': diff7_info_lower, 
                                'diff7_upper': diff7_info_upper, 
                                'pred_7_upper': pred_7_upper, 
                                'pred_7_lower': pred_7_lower, 
                                'pred7_med': pred_7_med,'as_extreme7': as_extreme7, 
                                'act7_upper': act_7_upper, 'act_7_lower': act_7_lower, 'act_7_med': act_7_med,
                                'early_sel_upper': early_sel_upper, 'early_sel_lower': early_sel_lower, 'early_sel_med': early_sel_med,
                                'late_sel_upper': late_sel_upper, 'late_sel_lower': late_sel_lower, 'late_sel_med': late_sel_med,

                                })
    df['sample_early1'] = early1
    df['sample_early2'] = early2
    df['sample_late1'] = late1
    df['sample_late2'] = late2

    return df
       



def get_main(metadata,species_dir,species, parent_samples, child_samples, freq_filtered, depth_filtered,early_interval = (1,2), late_interval = (6,7)):
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

    parent1_info = get_bootstrap_sel_coeffs(metadata,freq_children, med_depth_children, depth_children, reads_children, reads_children_opp, parent1_snps,parent2_snps, inoculumn, n_bootstraps = 1000)
  #  parent2_info = get_bootstrap_sel_coeffs(metadata,freq_children, med_depth_children, reads_children, reads_children_opp, parent2_snps, n_bootstraps = 1000)

    return pd.concat([count_parent1, count_parent2],axis=1), parent1_info #parent2_info 


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
    metadata = pd.read_csv('workflow/analysis/e003_coalescence_metadata_round4.csv')
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
    final_sample = 2
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
        early = (1,2)
        late=(6,7)
        _, parent1_info = get_main(metadata,species_dir,args.species, parent_samples, child_samples, 
            freq_filtered_in[parent_samples+child_samples],  depth_filtered_in[parent_samples+child_samples],early_interval = early, late_interval = late)

        inoculumn = ''.join(inoculumn.split('/'))
    
      #  count_parents.to_csv(f'{save_dir}/{inoculumn}_count_parents.csv')
        parent1_info.to_csv(f'{save_dir}/{inoculumn}_{str(early[0])}_{str(early[1])}_{str(late[0])}_{str(late[1])}_parent1_info_shift.csv')
       # parent2_info.to_csv(f'{save_dir}/{inoculumn}_{str(sample_init)}_{str(final_sample)}_parent2_info_shift.csv')

      #  distinguishing_snps.to_csv(f'{save_dir}/{inoculumn}_distinguishing_snps.csv')
       # freq_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_freq.csv.gz',compression = 'gzip')
        #depth_filtered_in.loc[distinguishing_snps.index.values,:].to_csv(f'{species_dir}/{inoculumn}_distinguishing_snps_depth.csv.gz',compression='gzip')
    


       
        
    

