from os.path import join
import pandas as pd
import os
from glob import glob
import numpy as np
configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_csv('config/sample_fnamesr3.csv')
df3=pd.read_csv('config/sample_fnames_round2.csv')
df2 = pd.read_csv('sample_fnames.csv')
df = pd.read_csv('sample_fnamesr4.csv')
# Parse sample names from df and generate sample list 
#df['SampleLane'] = df['Sample'].transform(lambda x: f'{x}_L002')
samples=list(df['SampleLane'].values.astype(str))
samples2=list(df2['Sample'].values.astype(str))
samples3 = list(df3['Sample'].values.astype(str))
samples = samples + samples2 + samples3
samples = list(df['SampleLane'].values.astype(str))
#print(samples)
df = df.set_index('SampleLane')
samples.remove('B_G2_G3_ACPP_AE_mBHI_mGAM_2_S170')
samples.remove('D_A5_A8_AA_AA_mGAM_mBHI_6_S293')
samples.remove('A_H7_D10_AE_AF_mGAM_mGAM_2_S91')
samples.remove('C_A12_E5_AF_AA_mBHI_mBHI_4_S204')
samples.remove('C_B12_e003Assembly_4_S216')
samples.remove('B_C2_G3_ACPP_AE_mBHI_mBHI_2_S122')
#print(df.head())
#samples = ['A5-e003Coalescence-mBHI-p7']
#f = df.loc[df['SampleLane'] != 'D4-e003Coalescence-mBHI-p5_S175_L003',:]
#df = df.loc[df['Sample'].isin(samples),:]

#samplelanes = list(set(df['SampleLane'].tolist()))
#df = df.set_index('SampleLane')
#print(samplelanes)
#print(samples)
#print(df)
#print(df.index.values)
#df = pd.read_csv('workflow/out/midas2_output/species/species_prevalence.tsv',delimiter= '\t').sort_values(by='sample_counts',ascending = False)
#good_species = df['species_id'].to_list()
#fnames1 = glob('workflow/out/midas2_output/merge/snps/*/*.snps_freqs.tsv')
#fnames2 = glob('workflow/out/midas2_output/merge/snps/*/*.snps_freqs.tsv.lz4')
#good_species = [f.split('/')[-1].split('.')[0] for f in fnames2] #+ [f.split('/')[-1].split('.')[0] for f in fnames2] 
#print(good_species)

def get_species_list():
    df = pd.read_csv('workflow/out/midas2_output/species/species_marker_median_coverage.tsv', delimiter = '\t').set_index('species_id')
    dfgood = df>= 5
    #dfgood = dfgood.sum(axis = 1)
    print(dfgood)    
    dfunique = pd.read_csv('workflow/out/midas2_output/species/species_unique_fraction_covered.tsv', delimiter = '\t').set_index('species_id') 
    dfgoodunique = dfunique>= .5
    dfgood = dfgood*dfgoodunique 
    dfgood = dfgood.sum(axis=1)
 #   fnames = glob('workflow/out/midas2_output/B4-e004Assembly-mBHI-p3/snps/*.snps.tsv.lz4')
  #  species = [fname.split('/')[-1].split('.')[0] for fname in fnames]
   # return species
    return list(dfgood[dfgood>5].index.values)

#print(samples[:10])
#print(df.index.values[:10])
#print('G4-e003Coalescence-mBHI-p5_S27' in samples)
species_list = [101346,102478, 100099]
fnames = glob('workflow/out/midas2_output/merge_10*/snps/10*/10*.snps_freqs.tsv*')
species_list = [fname.split('/')[-2] for fname in fnames]
#species_list.remove('102478')
#species_list.remove('101346')
#species_list.remove('100099')
species_list.remove('103686') #actuallly remove
#print('G4-e003Coalescence-mBHI-p5_S27' in df.index.values)
#species_list2 = [103117,100013, 101349,  101346,102478, 100099, 102279, 100146, 102438, 100196, 102506, 100120, 102528, 101294, 102327]

#species_list2 = [102506,102528,101349,  100196,102327]
#species_list = [100099, 101346,102320]
species_list2 = []
for species in species_list2: 
    print(species)
    species_list.remove(str(species))
#species_list.remove('103439') # add back
 #      'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
  #     'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
   #    'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
    #   'AE-AF-mGAM', 'AE-AF-mBHI',]
#samples.remove('H5-e003Coalescence-mGAM-p7_S184') # ADD BACK AFTER FIRST QUICK CHECK 
#samples.remove('Coalescence-F4-E9-AA-ACPP-mGAM-mBHI-3_S683')
#samples.remove('Coalescence-C4-A3-AA-ACPP-mBHI-mBHI-5_S647')
#species_list = species_list2
#print(df.head())
#print('wee',np.sort(samples))
#species_list = [101346, 100196]#,102478, 100099]
#species_list = [101346,100196,102478,100099, 102506, 100120]
#species_list = [101346]
#species_list.remove('101059')
#species_list.remove('102279')
#species_list.remove('103656')

#species_list = [100099,101059, 102279, 103656]
#species_list = [101346, 100196,  102506, ]
species_list = ['100146',
 '100910',
 '102478',
 '102438',
 '100196',
 '103117',
 '100074',
 '103439',
 '101958',
 '100044',
 '102528',
 '101349',
 '102506',
 '101346',
 '100099',
 '100120',
 '101400',
 '101294',
 '102057',
 '102945',
 '102327',
 '100013',
 '102320']
#print(species_list)
#species_list = [100146,102478,102478,102478,102478,102478,102478,102478,102478,102478,102478,100196, 100099]
#species_list = [100146,102478,100196,102438,100099,101346,102506]
print(len(samples))
#species_list = ['100196','100146']
#species_list = get_species_list()
#samples=['D_A4_A5_AF_AA_mBHI_mBHI_6_S292']
print(len(species_list))
#species_list = ['100099', '100146', '100196', '100910', '101059', '101294',
 #      '101346', '102438', '102478', '102506', '102528', '102544',
  #     '103117', '103656']
#<<<<<<< HEAD
species_list=[100196]
species_list=[100002, 100003, 100013, 100022, 100023, 100031, 100035, 100038,
       100044, 100059, 100071, 100074, 100078, 100084, 100087, 100088,
       100113, 100130, 100142, 100143, 100150, 100162, 100174, 100177,
       100181, 100187, 100190, 100193, 100198, 100199, 100205, 100212,
       100214, 100220, 100229, 100233, 100243, 100251, 100255, 100258,
       100262, 100293, 100296, 100646, 101029, 101059, 101060, 101141,
       101229, 101300, 101310, 101313, 101315, 101337, 101338, 101342,
       101345, 101361, 101362, 101367, 101378, 101382, 101400, 101433,
       101444, 101454, 101456, 101493, 101500, 101562, 101598, 101611,
       101638, 101694, 101714, 101747, 101958, 102040, 102143, 102272,
       102279, 102286, 102295, 102298, 102302, 102316, 102321, 102340,
       102347, 102353, 102445, 102448, 102455, 102470, 102492, 102528,
       102538, 102544, 102545, 102549, 102570, 102619, 102625, 102640,
       102684, 102720, 102794, 102876, 102896, 102903, 102926, 102945,
       103141, 103166, 103190, 103656, 103681, 103683, 103686, 103695,
       103727, 103772, 103876, 103886, 103937, 104099, 104552]
species_list = [100002, 100003, 100013, 100022, 100023, 100031, 100035, 100038,
       100044, 100057, 100059, 100071, 100074, 100084, 100087, 100088,
       100113, 100130, 100142, 100143, 100150, 100162, 100174, 100177,
       100181, 100187, 100193, 100198, 100199, 100205, 100212, 100214,
       100220, 100229, 100233, 100243, 100251, 100255, 100258, 100262,
       100293, 100296, 100646, 100758, 100818, 101029, 101059, 101060,
       101141, 101229, 101274, 101288, 101310, 101313, 101315, 101337,
       101338, 101342, 101345, 101361, 101362, 101367, 101378, 101382,
       101400, 101433, 101444, 101454, 101456, 101493, 101500, 101562,
       101598, 101611, 101638, 101694, 101714, 101747, 101958, 102035,
       102040, 102143, 102279, 102281, 102286, 102295, 102298, 102302,
       102316, 102321, 102323, 102340, 102347, 102353, 102445, 102448,
       102455, 102470, 102492, 102515, 102528, 102538, 102544, 102545,
       102549, 102570, 102619, 102625, 102640, 102645, 102684, 102794,
       102876, 102903, 102926, 102945, 103141, 103188, 103190, 103656,
       103681, 103683, 103686, 103695, 103727, 103772, 103796, 103876,
       103886, 103902, 103937, 104099, 104321, 104552]
#species_list.remove(100078)
#species_list.remove(100190)
species_list.remove(100057)

species_list.remove(100758)
species_list = [100002, 100003, 100013, 100022, 100023, 100031, 100035, 100038,
       100044, 100057, 100059, 100071, 100074, 100084, 100087, 100088,
       100099, 100113, 100120, 100130, 100142, 100143, 100146, 100150,
       100162, 100174, 100177, 100181, 100187, 100193, 100196, 100198,
       100199, 100205, 100212, 100214, 100220, 100229, 100233, 100243,
       100251, 100255, 100258, 100262, 100293, 100296, 100646, 100910,
       101029, 101059, 101060, 101141, 101229, 101288, 101294, 101310,
       101313, 101315, 101337, 101338, 101342, 101345, 101346, 101349,
       101361, 101362, 101367, 101378, 101382, 101400, 101433, 101444,
       101454, 101456, 101493, 101500, 101562, 101598, 101611, 101638,
       101694, 101714, 101747, 101958, 102040, 102057, 102143, 102279,
       102281, 102286, 102295, 102298, 102302, 102316, 102320, 102321,
       102327, 102340, 102347, 102353, 102438, 102445, 102448, 102455,
       102470, 102478, 102492, 102506, 102528, 102538, 102544, 102545,
       102549, 102570, 102619, 102625, 102640, 102645, 102684, 102794,
       102876, 102903, 102926, 102945, 103117, 103141, 103188, 103190,
       103439, 103656, 103681, 103683, 103686, 103695, 103727, 103772,
       103796, 103876, 103886, 103937, 104099, 104552]
species_list.remove(100057)
species_list.remove(101288)
species_list = [100002, 100003, 100013, 100022, 100023, 100031, 100035, 100038,
       100044, 100059, 100071, 100074, 100084, 100087, 100088, 100099,
       100113, 100120, 100130, 100142, 100143, 100146, 100150, 100162,
       100174, 100177, 100181, 100187, 100193, 100196, 100198, 100199,
       100205, 100212, 100214, 100220, 100229, 100233, 100243, 100251,
       100255, 100258, 100262, 100293, 100296, 100646, 100910, 101029,
       101059, 101060, 101141, 101229, 101294, 101310, 101313, 101315,
       101337, 101338, 101342, 101345, 101346, 101349, 101361, 101362,
       101367, 101378, 101382, 101400, 101433, 101444, 101454, 101456,
       101493, 101500, 101562, 101598, 101611, 101638, 101694, 101714,
       101747, 101958, 102040, 102057, 102143, 102279, 102286, 102295,
       102298, 102302, 102316, 102320, 102321, 102327, 102340, 102347,
       102353, 102438, 102445, 102448, 102455, 102470, 102478, 102492,
       102506, 102528, 102538, 102544, 102545, 102549, 102570, 102619,
       102625, 102640, 102684, 102794, 102876, 102903, 102926, 102945,
       103117, 103141, 103190, 103439, 103656, 103681, 103683, 103686,
       103695, 103727, 103772, 103796, 103876, 103886, 103937, 104099,
       104552]
#=======
species_list = ['103656']
#>>>>>>> 048bc47ef37a297b645a22b183fe842bd2d8d4a3
rule all:
    input:
         #expand("workflow/out/trimmed/{sample}-trimmed-pair1.fastq.gz",sample=samples),
        # expand("workflow/out/filter/{sample}-filtered.1.fastq.gz",sample=samples),
   #     expand("workflow/out/concat/{sample}-filtered.1.fastq.gz",sample=samples),
         #expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
         #"workflow/out/midas2_output/mergev4/species/species_prevalence.tsv",
       #  expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
        #"workflow/out/midas2_output/merge/snps/snps_summary.tsv",
      # "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
        # expand("workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv", species=species_list),
      #  expand("workflow/out/midas2_output/merge/snps/{species}/{species}.snps_freqs.tsv.gz", species=species_list),
        # expand("workflow/report/calculateDiversityDepthv3/{species}/{species}_diversity_df1.csv",species=species_list), 
        # expand("workflow/report/getSharedAlts/{species}/{species}_shared_80.csv",species=species_list),
         #expand("workflow/report/getSharedAlts/{species}/{species}_shared_20.csv",species=species_list)    
#=======
      #   expand("workflow/report/getSharedAlts/{species}/{species}_shared_80.csv")
       # expand("workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv",species=species_list),  
     #   expand("workflow/report/track_snpsv2_ALL/{species}/done.txt", species=species_list),   
      #   expand("workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/done.txt", species=species_list), 
 #        expand("workflow/report/calculateFixedDiffsFastv3/{species}/{species}_fixed_diffs.csv",species=species_list),
         expand("workflow/report/track_snpsv3_sel_bootstrapv3/{species}/done.txt",  species=species_list),
         expand( "workflow/report/track_snpsv5_shift_self_test_mod_shift/{species}/done12.txt",species=species_list),
#"workflow/report/track_snpsv3_sel_bootstrapv3/{species}/done.txt"
        # expand("workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done.txt",  species=species_list),
        # expand("workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done23.txt",  species=species_list)
       #  expand("workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done01.txt",  species=species_list),
        # expand("workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done12.txt",  species=species_list),
         #expand("workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done23.txt",  species=species_list),
         #expand("workflow/report/track_snpsv2_ALL_bootstrapv3_both/{species}/done.txt",species=species_list)
       # expand("workflow/report/track_snpsv2_ALL_same_subject/{species}/done.txt",species=species_list),
      #  workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv
#  "workflow/report/calculateFixedDiffs/100013/100013_fixed_diffs.csv"
 # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),

#include: "workflow/rules/processRawReads_no_concatenation.smk",
#include: "workflow/rules/processRawReads.smk",
#include: "workflow/rules/runMIDAS2.smk",
#include: "workflow/rules/runMIDAS2_population.smk"
include: "workflow/rules/processMIDAS2.smk"
#includ: "workflow/rules/processSNPs.smk"
