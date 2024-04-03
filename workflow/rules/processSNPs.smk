import pandas as pd 

# Run the MIDAS species module to generate species profiles.



rule compute_populationSNVs:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        "workflow/out/midas2_output/merge/snps/snps_summary.tsv"
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"],
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_snps --samples_list workflow/out/list_of_samples.tsv  --midasdb_name {params.midasdb} \
            --advanced --midasdb_dir {params.midasdb_dir} --genome_coverage 0.7 --num_cores {threads} workflow/out/midas2_output/merge
        """

rule get_filtered_dfs:
    input:
        snpsDepth="/workflow/out/midas2_output/merge/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="/workflow/out/midas2_output/merge/snps/{species}/{species}.snps_freq.tsv",
        snpsInfo="/workflow/out/midas2_output/merge/snps/{species}/{species}.snps_info.tsv",
        snpsSummary=config['MIDASRundir']+"/workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_summary.txt.gz"
    output:
        "workflow/out/filtered_snps/{species}/diversity_df.csv",
        "workflow/out/filtered_snps/{species}/median_depth_df.csv"
    params:
        outdir="workflow/out/midas2_output/merge/snps/",
        savedir="workflow/out/filtered_snps/",
        species="{species}",
        
    conda:
        "../../workflow/envs/python3-analysis.yml"
    shell:
        """
        python3 workflow/scripts/analyzeMIDAS/get_filtered_subject_dfs.py {params.outdir} {params.savedir} {params.species} 
        """

