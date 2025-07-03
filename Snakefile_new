from os.path import join
import pandas as pd
import os
from glob import glob
import numpy as np
configfile: "config/config.yaml"

# Convert list of samples to a dataframe.
df=pd.read_csv('config/sample_fnamesr3.csv')
samples = df['SampleLane'].values #+ samples2
print(np.sort(samples))
df = df.set_index('SampleLane')
print('B10-e003Coalescence-mGAM-p1_S769' in samples)
rule all:
    input:
        expand("workflow/out/trimmed/{sample}-trimmed-pair1.fastq.gz",sample=samples),
        expand("workflow/out/filter/{sample}-filtered.1.fastq.gz",sample=samples),
   #     expand("workflow/out/concat/{sample}-filtered.1.fastq.gz",sample=samples),
        expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples),
      # "workflow/out/midas2_output/mergev2/species/species_prevalence.tsv",
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples),
#        "workflow/out/midas2_output/merge/snps/snps_summary.tsv",
      # "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
      #  expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples),
       # expand("workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz", species=species_list),
      #  expand("workflow/out/midas2_output/merge/snps/{species}/{species}.snps_freqs.tsv.gz", species=species_list),
      #  expand("workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df1.csv",species=species_list),     
       # expand("workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv",species=species_list),  
       # expand("workflow/report/track_snpsv2/{species}/done.txt", species=species_list),   
       # expand("workflow/report/track_snpsavg/{species}/done.txt", species=species_list), 
        #expand("workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv",species=species_list),
      #  workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv
#  "workflow/report/calculateFixedDiffs/100013/100013_fixed_diffs.csv"
 # "workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        #"workflow/out/midasOutput/species/abundantSpecies.txt",
       # expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt", subject=subjects),

include: "workflow/rules/processRawReads_no_concatenation.smk",
#include: "workflow/rules/processRawReads.smk",
include: "workflow/rules/runMIDAS2.smk",
#include: "workflow/rules/runMIDAS2_population.smk"
#include: "workflow/rules/processMIDAS2.smk"
#include: "workflow/rules/processSNPs.smk"
