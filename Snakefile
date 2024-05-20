from os.path import join
import pandas as pd
import os

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

rule all:
    input:
 #       expand("workflow/out/trimmed/{samplelane}-trimmed-pair1.fastq.gz",samplelane=samplelanes),
  #      expand("workflow/out/filter/{samplelane}-filtered.1.fastq.gz",samplelane=samplelanes),
   #     expand("workflow/out/concat/{sample}-filtered.1.fastq.gz",sample=samples),
        expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
        "workflow/out/midas2_output/species/abundant_species.csv",
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
        "workflow/out/midas2_output/merge/snps/snps_summary.tsv",
       # "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
       # "workflow/out/midas2_output/merge/snps/100013/100013.snps_freqs.tsv.gzip",
        "workflow/report/calculateDiversityDepth/100013/100013_diversity_df.csv",      
 # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),


#include: "workflow/rules/processRawReads.smk",
include: "workflow/rules/runMIDAS2.smk",
include: "workflow/rules/processMIDAS2.smk"
#include: "workflow/rules/processSNPs.smk"
