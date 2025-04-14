# Analyze the SNP calls for each species.
# Process the SNP calls, plot the SFS,
# and calculate the number of fixed differences between pairs of samples.

rule compressgz: 
    input:
        snpsDepth="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_info.tsv",
    output:
        snpsDepth="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
        snpsInfo="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_info.tsv.gz",
    shell:
        """
        gzip {input.snpsDepth}
        gzip {input.snpsFreq}
        gzip {input.snpsInfo}
        """

rule calculateDiversity:
    input:
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateDiversityDepthv3/{species}/{species}_diversity_df1.csv" 
    params:
        indir="workflow/out/midas2_output/mergev3_{species}/snps/",
        outdir="workflow/report/calculateDiversityDepthv3/",
        species="{species}"
#    conda:
 #       "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_diversity_med_depth_basic_filtering.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"


rule calculateFixedDiffs:
    input:
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateFixedDiffs/{species}/{species}_fixed_diffs.csv"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/calculateFixedDiffs/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_pairwise_fixed_diffs_only_in.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"

rule calculateFixedDiffsFast:
    input:
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
    output:
        "workflow/report/calculateFixedDiffsFastv3/{species}/{species}_fixed_diffs.csv"
    params:
        indir="workflow/out/midas2_output/mergev3_{species}/snps/",
        outdir="workflow/report/calculateFixedDiffsFastv3/",
        species="{species}"
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        "python3 workflow/scripts/get_pairwise_fixed_diffs_only_fast.py --outdir {params.outdir} --indir {params.indir} --species {params.species}"



rule trackSNPs:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/track_snpsv2/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """

rule trackSNPsAVG:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_ALL/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/track_snpsv2_ALL/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


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



rule trackSNPsAVG_bootstrap:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_ALL_bootstrap/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/track_snpsv2_ALL_bootstrap/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_bootstrap.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


rule trackSNPsAVG_bootstrapv3:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev3_{species}/snps/",
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
        snpsFreq="workflow/out/midas2_output/mergev2_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_sel_bootstrap/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev2_{species}/snps/",
        outdir="workflow/report/track_snpsv2_sel_bootstrap/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_bootstrap_pairs.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """



rule trackSNPsAVG_shift_sel:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_shift_self_test_13/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev3_{species}/snps/",
        outdir="workflow/report/track_snpsv2_shift_self_test_13/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_sample_quads_pairs13.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """


rule trackSNPsAVG_shift_sel_mod:
    input:
       # snpsDepth="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_depth.tsv.gz",
        snpsFreq="workflow/out/midas2_output/mergev3_{species}/snps/{species}/{species}.snps_freqs.tsv.gz",
       # snpsInfo="workflow/out/midas2_output/merge_{species}/snps/{species}/{species}.snps_info.tsv.gz",
       # wo="workflow/report/calculateDiversityDepth/{species}/{species}_diversity_df.csv"
    output:
        "workflow/report/track_snpsv2_shift_self_test_mod_13/{species}/done.txt"
    params:
        indir="workflow/out/midas2_output/mergev3_{species}/snps/",
        outdir="workflow/report/track_snpsv2_shift_self_test_mod_13/",
        #species={species}
  #  conda:
   #     "../../workflow/envs/snps_analysis_tools-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/track_snps_avg_v2_sample_quads_pairs_modular.py --outdir {params.outdir} --indir {params.indir} --species {wildcards.species}
        touch {params.outdir}/{wildcards.species}/done.txt
        """
