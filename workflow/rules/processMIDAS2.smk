# Analyze the SNP calls for each species.
# Process the SNP calls, plot the SFS,
# and calculate the number of fixed differences between pairs of samples.

rule decompresslz4: 
    input:
        snpsDepth="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_depth.tsv.lz4",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.lz4",
        snpsInfo="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_info.tsv.lz4",
    output:
        snpsDepth="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_info.tsv",
    shell:
        """
        lz4 -d --rm {input.snpsDepth} {output.snpsDepth}
        lz4 -d --rm {input.snpsFreq} {output.snpsFreq}
        lz4 -d --rm {input.snpsInfo} {output.snpsInfo}
        """

rule decompresslz4: 
    input:
        snpsDepth="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_info.tsv",
    output:
        snpsDepth="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.ts.gz",
        snpsInfo="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_info.tsv.gz",
    shell:
        """
        gzip {input.snpsDepth} {output.snpsDepth}
        gzip {input.snpsFreq} {output.snpsFreq}
        gzip {input.snpsInfo} {output.snpsInfo}
        """

rule calculateDiversity:
    input:
        snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
        snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
    output:
        "workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    params:
        indir="workflow/out/midas2_output/merge_{species}/snps/",
        outdir="workflow/report/calculateDiversityDepth/",
        species="{species}"
#    conda:
 #       "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_diversity_med_depth_basic_filtering.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"


rule calculateFixedDiffs:
    input:
        snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
        snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
        wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv"
    params:
        indir="workflow/out/midas2_output/merge_{species}/snps/",
        outdir="workflow/report/calculateFixedDiffs/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_pairwise_fixed_diffs.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"




rule trackSNPs:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snps/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/merge_{species}/snps/",
        outdir="workflow/report/track_snps/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """

