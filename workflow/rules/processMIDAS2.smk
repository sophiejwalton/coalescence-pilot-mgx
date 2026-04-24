# Analyze the SNP calls for each species.
# Process the SNP calls, plot the SFS,
# and calculate the number of fixed differences between pairs of samples.

rule decompresslz4:
    input:
        snpsDepth="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv.lz4",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.lz4",
        snpsInfo="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv.lz4",
    output:
        snpsDepth="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv",
    shell:
        """
        lz4 -d {input.snpsDepth} {output.snpsDepth}
        lz4 -d {input.snpsFreq} {output.snpsFreq}
        lz4 -d {input.snpsInfo} {output.snpsInfo}
        """
rule compressgz: 
    input:
        snpsDepth="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv",
    output:
        snpsDepth="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
        snpsInfo="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv.gz",
    shell:
        """
        gzip {input.snpsDepth}
        gzip {input.snpsFreq}
        gzip {input.snpsInfo}
        """

rule calculateDiversity:
    input:
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateDiversityDepthv3/{species}/{species}_diversity_df3.csv" 
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/calculateDiversityDepthv3/",
        species="{species}"
#    conda:
 #       "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_diversity_med_depth_basic_filtering.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"


rule calculateFixedDiffs:
    input:
        snpsFreq="workflow/out/midas2v3_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv"
    params:
        indir="workflow/out/midas2v3_output/mergev2_{species}/snps/",
        outdir="workflow/report/calculateFixedDiffs/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_pairwise_fixed_diffs_only_in.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"

rule calculateFixedDiffsFast:
    input:
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateFixedDiffsFastv3/{species}/{species}_fixed_diffs.csv"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/calculateFixedDiffsFastv3/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_pairwise_fixed_diffs_only_fast.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"


rule getsharedalts:
    input:
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/getSharedAlts/{species}/{species}_shared_80.csv"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/getSharedAlts/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_shared_alts_only_fast.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"



rule trackSNPsAVG_same_subject:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_ALL_same_subject/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/track_snpsv2_ALL_same_subject/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_same_subject.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """



rule trackSNPsAVG_bootstrapv3:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/done.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_ALL_bootstrapv3/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_bootstrap.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


rule trackSNPsAVG_bootstrap_pairs:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_sel_bootstrapv3/{species}/done.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_sel_bootstrapv3/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_bootstrap_pairs.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


rule trackSNPsAVG_bootstrap_pairsv3:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv3_sel_bootstrapv3/{species}/done.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv3_sel_bootstrapv3/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v3_bootstrap_pairs.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


rule trackSNPsAVG_bootstrapv3both:
    input:
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/track_snpsv2_ALL_bootstrapv3_both/{species}/done.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_ALL_bootstrapv3_both/",
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_bootstrap_both.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


#workflow/report/track_snpsv2_shift_self_test_12_shift
rule trackSNPsAVG_shift_sel_mod_shift01:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done01.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_shift_self_test_mod_shift/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_sample_quads_pairs_modular_shift.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species} --interval 01
        touch {params.outdir}/{wildcards.species}/done01.txt
        """

#workflow/report/track_snpsv2_shift_self_test_12_shift
rule trackSNPsAVG_shift_sel_mod_shift12:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done12.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_shift_self_test_mod_shift/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_sample_quads_pairs_modular_shift.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species} --interval 12
        touch {params.outdir}/{wildcards.species}/done12.txt
        """


#workflow/report/track_snpsv2_shift_self_test_12_shift
rule trackSNPsAVG_shift_sel_mod_shift23:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2v3_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_shift_self_test_mod_shift/{species}/done23.txt"
    params:
        indir="workflow/out/midas2v3_output/mergevfinal_{species}/snps/",
        outdir="workflow/report/track_snpsv2_shift_self_test_mod_shift/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_sample_quads_pairs_modular_shift.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species} --interval 23
        touch {params.outdir}/{wildcards.species}/done23.txt
        """