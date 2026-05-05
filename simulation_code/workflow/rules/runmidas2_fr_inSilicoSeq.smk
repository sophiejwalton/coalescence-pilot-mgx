import pandas as pd 

def get_abundance_species():
    # pass
    df= pd.read_csv('/oak/stanford/groups/dpetrov/swalton/coalescence-pilot-mgx/workflow/out/midas2_output/species/abundant_species.csv')
    return ','.join(list(df['species_id'].astype(str)))

rule identifySNVs:
    input:
        r1='workflow/out/insilicoSeqOutput20ones/{trial_part}/insilico-reads.fastq_R1.fastq',
        r2='workflow/out/insilicoSeqOutput20ones/{trial_part}/insilico-reads.fastq_R2.fastq',
    output:
        profile="workflow/out/midasoutput_with_ones/{trial_part}/snps/snps_summary.tsv"
    params:
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir'],
        species_list=get_abundance_species()
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 run_snps --sample_name {wildcards.trial_part} -1 {input.r1} -2 {input.r2}  --midasdb_name {params.midasdb} \
            --midasdb_dir {params.midasdb_dir} --num_cores {threads} workflow/out/midasoutput_with_ones \
            --advanced --site_depth 1 --snp_maf 0.001 \
            --species_list {params.species_list} --select_threshold=-1
        """

rule compute_populationSNVs:
    input:
        expand("workflow/out/midasoutput_with_ones/{trial_part}/snps/snps_summary.tsv",trial_part=trial_parts)
    output:
        #"workflow/out/midas2_output/merge_v2/snps/snps_summary.tsv",
        "workflow/out/midasoutput_with_ones/mergeonly_two_genomes_100196/snps/100196/100196.snps_freqs.tsv.lz4",
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"],
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir'],
      #  species_list = get_species_list()
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_snps --samples_list 'workflow/out/list_of_samples.tsv' --species_list 100196  --midasdb_name {params.midasdb} --genome_depth 1 --site_depth 1  --site_prev 0.0 --snp_maf 0.01  --advanced --midasdb_dir {params.midasdb_dir} --snp_type any --genome_coverage 0.5 --num_cores {threads} workflow/out/midasoutput_with_ones/mergeonly_two_genomes_100196
        """
