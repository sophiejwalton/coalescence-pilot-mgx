import pandas as pd 

# Run the MIDAS species module to generate species profiles.
rule profileSpeciesAbundances:
    input:
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
    output:
        profile="workflow/out/midasOutput/{sample}/species/species_profile.txt"
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/MIDASpython2-no-builds.yml"
    shell:
        """
        run_midas.py species workflow/out/midasOutput/{wildcards.sample} \
            -1 {input.r1} -2 {input.r2} -t {threads} -d /home/groups/bhgood/midas_db_ONLYHUMAN/midas_db_v1.2
        """

rule get_combined_species:
    input:
        expand("workflow/out/midasOutput/{sample}/species/species_profile.txt",sample=samples)
    output:
        "workflow/out/midasOutput/species/species_profile_all_abundant.csv"
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"]
    shell:
        """
        python3 workflow/scripts/combine_species_profiles.py {params.minCoverage}
        """

# Helper function to extract the list of abundant species in a household
# based on the sample provided.
def getAbundantSpecies():
    df=pd.read_csv("workflow/out/midasOutput/species/species_profile_all_abundant.csv")
    species_list=list(df['species_id'].unique())
    species_list=",".join(species_list)
    return species_list


# Run the SNPs module on each sample using the set of species
# found to be abundant in the household that the sample belongs to.
rule callSNPs:
    input:
        splist="workflow/out/midasOutput/species/species_profile_all_abundant.csv",
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
    params:
        species_list = getAbundantSpecies()
    output:
        "workflow/out/midasOutput/{sample}/snps/summary.txt"
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/MIDASpython2-no-builds.yml"
    shell:
        """
        run_midas.py snps workflow/out/midasOutput/{wildcards.sample} \
            -1 {input.r1} -2 {input.r2} -t {threads} \
            --species_id {params}
        """