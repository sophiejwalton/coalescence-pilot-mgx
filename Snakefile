from os.path import join
import pandas as pd
import os
from glob import glob
configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_csv('sample_fnames.csv')

# Parse sample names from df and generate sample list 
#df['SampleLane'] = df['Sample'].transform(lambda x: f'{x}_L002')
samples=list(set(df['Sample'].sort_values().tolist()))
#samples = ['A5-e003Coalescence-mBHI-p7']
df = df.loc[df['SampleLane'] != 'D4-e003Coalescence-mBHI-p5_S175_L003',:]
#df = df.loc[df['Sample'].isin(samples),:]

samplelanes = list(set(df['SampleLane'].tolist()))
df = df.set_index('SampleLane')
#print(samplelanes)
#print(samples)
#print(df)
#print(df.index.values)
df = pd.read_csv('workflow/out/midas2_output/species/species_prevalence.tsv',delimiter= '\t').sort_values(by='sample_counts',ascending = False)
good_species = df['species_id'].to_list()
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

species_list = get_species_list()
#print(species_list)
#print(species_list)
#species_list.remove('100214')
#all_done = glob('workflow/out/midas2_output/merge/snps/*/*.snps_freqs.tsv.lz4')
#all_done = glob('workflow/out/midas2_output/merge_*/snps/*/*.snps_freqs.tsv')
#print(all_done)
print(species_list)
#for fname in all_done:
 #   species = fname.split('/')[-2]
  #  if species in species_list:
   # 	species_list.remove(fname.split('/')[-2]) 
#species_list.remove('101367')
#species_list.remove('101400') # make sure to add this one back in
#species_list.remove('103739')
#species_list.remove('100818')
#species_list.remove('102843')
#species_list.remove('102482')
#species_list.remove('101787')
#species_list.remove('101274')
#species_list.remove(102554)
species_list.remove(102544)
#species_list.remove('100672')
species_list.remove(100758)
species_list.remove(101315)
species_list.remove(101611)
species_list.remove(100144)
species_list.remove(102448) # add back in
species_list.remove(100002) # add back in 
species_list.remove(103686)
species_list.remove(100212)
species_list.remove(100193)
species_list.remove(101302) # add back in 
species_list.remove(103891) # add back in 
species_list.remove(103899)
species_list.remove(101310)
species_list.remove(100262) # add back in 
species_list.remove(100060)
species_list.remove(100033) # add back in 
species_list.remove(100263) # add back in 
species_list.remove(100249) # add back in 
species_list.remove(101456) # add back in 
species_list.remove(100231)
#species_list.remove('103326')
#species_list.remove('101021')
#species_list.remove('102289')
#species_list.remove('102890')
#species_list.remove('100293')
#species_list.remove('103902')
#species_list.remove('100279')
#species_list.remove('103188')
#species_list.remove('104321')
#species_list.remove('102515')
#species_list.remove('100057')
#species_list.remove('100196') # ADD BACK IN ROUND 2`
#species_list.remove('102327') # make sure to add this one back
#species_list.remove('102478') # add back in
#species_list.remove('101346') # add back in 
#species_list.remove('100099') #add back in 
#species_list.remove('100196') #add back in 
#species_list.remove('102506') # add back in
#species_list.remove('101367') # add back in 
print('YAAAAAYY')
print(len(species_list))
#species_list.remove('101346')
#species_list.remove('102492')
#species_list.remove('103439')
#species_list.remove('100196')
#good_species.remove('101346')
inoculumns = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
rule all:
    input:
 #       expand("workflow/out/trimmed/{samplelane}-trimmed-pair1.fastq.gz",samplelane=samplelanes),
  #      expand("workflow/out/filter/{samplelane}-filtered.1.fastq.gz",samplelane=samplelanes),
   #     expand("workflow/out/concat/{sample}-filtered.1.fastq.gz",sample=samples),
   #     expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
  #      "workflow/out/midas2_output/species/abundant_species.csv",
       # expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
#        "workflow/out/midas2_output/merge/snps/snps_summary.tsv",
      # "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
 #       expand("workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_freqs.tsv", species=species_list),
      #  expand("workflow/out/midas2_output/merge/snps/{species}/{species}.snps_freqs.tsv.gz", species=species_list),
       # expand("workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv",species=species_list),     
       # expand("workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv",species=species_list),  
        expand("workflow/report/track_snps/{species}/done.txt", species=species_list)    
#  "workflow/report/calculateFixedDiffs/100013/100013_fixed_diffs.csv"
 # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),


#include: "workflow/rules/processRawReads.smk",
#include: "workflow/rules/runMIDAS2.smk",
include: "workflow/rules/processMIDAS2.smk"
#include: "workflow/rules/processSNPs.smk"
