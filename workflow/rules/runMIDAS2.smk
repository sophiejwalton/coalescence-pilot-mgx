import pandas as pd 

# Run the MIDAS species module to generate species profiles.
rule identifySpecies:
    input:
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
    output:
        profile="workflow/out/midas2_output/{sample}/species/species_profile.tsv"
    params:
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 run_species --sample_name {wildcards.sample} -1 {input.r1} -2 {input.r2} --midasdb_name {params.midasdb} \
            --midasdb_dir {params.midasdb_dir} --num_cores {threads} workflow/out/midas2_output
        """

rule get_abundant_species:
    input:
       expand("workflow/out/midas2_output/{sample}/species/species_profile.tsv",sample=samples)
    output:
        profile="workflow/out/midas2_output/species/abundant_species.csv"
    params:
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        python3 workflow/scripts/get_abundant_species.py 
        """

def get_abundance_species():
#     pass
    df= pd.read_csv('workflow/out/midas2_output/species/abundant_species.csv')
    return ','.join(list(df['species_id'].astype(str)))

rule identifySNVs:
    input:
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz"),
        good_species='workflow/out/midas2_output/species/abundant_species.csv',
        species="workflow/out/midas2_output/{sample}/species/species_profile.tsv"
    output:
        profile="workflow/out/midas2_output/{sample}/snps/snps_summary.tsv"
    params:
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir'],
        species_list=get_abundance_species()
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 run_snps --sample_name {wildcards.sample} -1 {input.r1} -2 {input.r2}  --midasdb_name {params.midasdb} \
            --midasdb_dir {params.midasdb_dir} --num_cores {threads} workflow/out/midas2_output  \
            --advanced \
            --species_list {params.species_list} --select_threshold=-1
        """

rule compute_populationSNVs:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        "workflow/out/midas2_output/merge/snps_summary.tsv"
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


rule compute_populationSNVs_prominent:
    input:
        expand("workflow/out/midas2_output/{sample}/snps/snps_summary.tsv",sample=samples)
    output:
        "workflow/out/midas2_output/merge_bacteroides/snps_summary.tsv"
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"],
        midasdb=config['midasdb'],
        midasdb_dir=config['midasdb_dir']
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/midas2_sw-no-builds.yml"
    shell:
        """
        midas2 merge_snps --samples_list workflow/out/list_of_samples.tsv  --midasdb_name {params.midasdb} --species_list 101346,102478 \
            --advanced --midasdb_dir {params.midasdb_dir} --genome_coverage 0.7 --num_cores {threads} workflow/out/midas2_output/merge_bacteroides
        """

