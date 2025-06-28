import pandas as pd 

rule mergeSpecies:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        "workflow/out/midas2_output/mergev2/species/species_prevalence.tsv"
    params:
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_species --samples_list workflow/out/list_of_samples.tsv --min_cov 1 midas2_output/mergev2
        """


rule compute_populationSNVs:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        #"workflow/out/midas2_output/merge_v2/snps/snps_summary.tsv",
        "workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.lz4",
        "workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv.lz4",
        "workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv.lz4",
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"],
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir'],
        species_list = get_species_list()
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_snps --samples_list workflow/out/list_of_samples.tsv --species_list {wildcards.species}  --midasdb_name {params.midasdb} --robust_chunk --genome_depth 1 --site_depth 1  --site_prev 0.0 --snp_maf 0.01  --advanced --midasdb_dir {params.midasdb_dir} --snp_type any --genome_coverage 0.5 --num_cores {threads} workflow/out/midas2_output/mergevfinal_{wildcards.species}
        """


rule compute_populationSNVs_prominent:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        "workflow/out/midas2_output/merge_bacteroides/snps/snps_summary.tsv",
        "workflow/out/midas2_output/merge_b/snps/101346/101346.snps_freqs.tsv.lz4" 
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"],
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_snps --samples_list workflow/out/list_of_samples.tsv  --midasdb_name {params.midasdb} --species_list 101346 --site_depth 1  --site_prev 0.0 --snp_maf 0.01  --advanced --midasdb_dir {params.midasdb_dir} --snp_type any --genome_depth 10 --genome_coverage 0.5  --robust_chunk  --num_cores {threads} workflow/out/midas2_output/merge_bacteroides
        """


rule decompresslz4: 
    input:
        snpsDepth="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv.lz4",
        snpsFreq="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv.lz4",
        snpsInfo="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv.lz4",
    output:
        snpsDepth="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_depth.tsv",
        snpsFreq="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_freqs.tsv",
        snpsInfo="workflow/out/midas2_output/mergevfinal_{species}/snps/{species}/{species}.snps_info.tsv",
    shell:
        """
        lz4 -d {input.snpsDepth} {output.snpsDepth}
        lz4 -d {input.snpsFreq} {output.snpsFreq}
        lz4 -d {input.snpsInfo} {output.snpsInfo}
        """
