from os.path import join
import pandas as pd
import os
from glob import glob
configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_csv('config/sample_fnames_round2.csv')

# Parse sample names from df and generate sample list 
#df['SampleLane'] = df['Sample'].transform(lambda x: f'{x}_L002')
samples=list(set(df['Sample'].sort_values().tolist()))
df = df.set_index('Sample')
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


inoculumns = ['AA-AE-mGAM', 'AA-AF-mGAM', 
       'AA-AC/PP-mGAM', 'AA-AC/PP-mBHI', 'AA-AE-mBHI', 'AA-AF-mBHI',
       'AC/PP-AE-mGAM', 'AC/PP-AF-mGAM', 
       'AC/PP-AE-mBHI', 'AC/PP-AF-mBHI', 
       'AE-AF-mGAM', 'AE-AF-mBHI',
     ]
rule all:
    input:
        expand("workflow/out/trimmed/{sample}-trimmed-pair1.fastq.gz",sample=samples),
        expand("workflow/out/filter/{sample}-filtered.1.fastq.gz",sample=samples),
   #     expand("workflow/out/concat/{sample}-filtered.1.fastq.gz",sample=samples),
   #     expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
  #      "workflow/out/midas2_output/species/abundant_species.csv",
       # expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
#        "workflow/out/midas2_output/merge/snps/snps_summary.tsv",
      # "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
     #   expand("workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_freqs.tsv", species=species_list),
      #  expand("workflow/out/midas2_output/merge/snps/{species}/{species}.snps_freqs.tsv.gz", species=species_list),
       # expand("workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv",species=species_list),     
       # expand("workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv",species=species_list),  
#        expand("workflow/report/track_snps/{species}/done.txt", species=species_list)    
#  "workflow/report/calculateFixedDiffs/100013/100013_fixed_diffs.csv"
 # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),


include: "workflow/rules/processRawReads_no_concatenation.smk",
#include: "workflow/rules/processRawReads.smk",
#include: "workflow/rules/runMIDAS2.smk",
#include: "workflow/rules/processMIDAS2.smk"
#include: "workflow/rules/processSNPs.smk"
