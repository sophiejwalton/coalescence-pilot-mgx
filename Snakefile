from os.path import join
import pandas as pd
import os

configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_csv('config/sample_fnames.csv')

# Parse sample names from df and generate sample list 
#df['SampleLane'] = df['Sample'].transform(lambda x: f'{x}_L002')
samples=list(set(df['Sample'].tolist()))
samplelanes = list(set(df['SampleLane'].tolist()))
df = df.set_index('Sample')
#print(df)
#print(samples)
#print(samplelanes)
#print(df['trim1'].values)
#print(df['trim2'].values)
#print(df['read1'].values)
#print(df['read2'].values)
#samples.remove('AD_0809_SW_26')
# Parse the list of species analyzed in the MIDAS snps module.
#if(os.path.isdir("workflow/out/midasOutput/snps")):
    # Iterate through the species analyzed by the snps module.
   # dirs=[x[0] for x in os.walk("workflow/out/midasOutput/snps/HouseholdTransmission-Stool")]
    # Remove the first element, which is the directory itself without subdirectories.
   # dirs.pop(0)
    # Parse the species names and generate a list.
  #  snpsSpecies=[]
   # for species in dirs:
      #  snpsSpecies.append(species.split("/")[5])

rule all:
    input:
        expand("workflow/out/trimmed/{samplelane}-trimmed-pair1.fastq.gz",samplelane=samplelanes),
        expand("workflow/out/filter/{sample}-filtered.1.fastq.gz",sample=samples),
        expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
        "workflow/out/midas2_output/species/abundant_species.csv",
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
        "workflow/out/midas2_output/merge/snps/snps_summary.tsv",
        "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
       # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),


include: "workflow/rules/processRawReads.smk",
include: "workflow/rules/runMIDAS2.smk",

