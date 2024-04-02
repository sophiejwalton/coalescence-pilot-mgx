# Run the MIDAS species module to generate species profiles.
rule profileSpeciesAbundances:
    input:
        r1=join(config["filterdir"],"{samplelane}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{samplelane}-filtered.2.fastq.gz")
    output:
        profile="workflow/out/midasOutput/{samplelane}/species/species_profile.txt"
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/MIDASpython2-no-builds.yml"
    shell:
        """
        run_midas.py species workflow/out/midasOutput/{wildcards.samplelane} \
            -1 {input.r1} -2 {input.r2} -t {threads}
        """

# Annotate each species profile with the sample name for further processing.
rule annotateSpeciesAbundancesBySubject:
    input:
        "workflow/out/midasOutput/{sample}/species/species_profile.txt"
    output:
        "workflow/out/midasOutput/{sample}/species/species_profile_subject.txt"
    shell:
        "sed 's/$/\t{wildcards.sample}/g' {input} > {output}"

# Concatenate all species abundance profiles into a single file.
rule concatenateSpeciesAbundances:
    input:
        expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/species/species_profile_subject.txt",sample=samples)
    output:
        "workflow/out/midasOutput/species/species_profile_all.txt"
    shell:
        # Concatenate all species profiles without the file headers.
        "cat <( tail -n +2 {input} ) | grep -v '==>' | sed '/^$/d' | awk '$3 > 0' |"
        # Add the file header back and include the new "sample" field.
        "sed '1 i\species_id\tcount_reads\tcoverage\trelative_abundance\tsample'> {output}"

# Identify the abundant species to analyze for SNPs
rule identifyAbundantSpecies:
    input:
        "workflow/out/midasOutput/species/species_profile_all.txt"
    params:
        minCoverage=config["runMIDAS_speciesMinCoverage"]
    output:
        "workflow/out/midasOutput/species/abundantSpecies.txt"
    shell:
        # Retain species with coverage greater than the minimum coverage threshold.
        "cat <( tail -n +2 {input} ) | grep -v '==>' | sed '/^$/d' | awk '$3 > {params.minCoverage}"
            "' | cut -f1 | sort | uniq > {output}"


# Helper function to extract the list of abundant species in a household
# based on the sample provided.
def getAbundantSpecies():
   
    with open("".join(["workflow/out/midasOutput/species/abundantSpecies.txt"])) as f:
        species=f.read().splitlines()
    # Convert the list of abundant species into a string.
    species=",".join(species)
    return species

# Run the SNPs module on each sample using the set of species
# found to be abundant in the household that the sample belongs to.
rule callSNPs:
    input:
        r1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
        r2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz")
    params:
        lambda wildcards: getAbundantSpeciesFromSample(wildcards.sample)
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

# Merge SNPs from all samples in all households.
rule mergeSNPsAllHouseholds:
    input:
        expand("workflow/out/midasOutput/HouseholdTransmission-Stool-{sample}/snps/summary.txt",sample=samples)
    params:
        sampledirs=",".join("workflow/out/midasOutput/HouseholdTransmission-Stool-" + s for s in samples),
        outdir="workflow/out/midasOutput/snps/HouseholdTransmission-Stool"
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/MIDASpython2-no-builds.yml"
    shell:
        """
        mkdir -p {params.outdir}
        merge_midas.py snps {params.outdir} -t list -i {params.sampledirs} \
            --sample_depth 5 --site_depth 3 --min_samples 1 --max_species 150 \
            --site_prev 0.0 --threads {threads}
        """
"""

# Compress the output files from the MIDAS merge SNPs function.
rule compressSNPoutput:
    input:
        snpsDepth="{dir}/snps_depth.txt",
        snpsFreq="{dir}/snps_freq.txt",
        snpsInfo="{dir}/snps_info.txt",
        snpsSummary="{dir}/snps_summary.txt"
    output:
        snpsDepth="{dir}/snps_depth.txt.gz",
        snpsFreq="{dir}/snps_freq.txt.gz",
        snpsInfo="{dir}/snps_info.txt.gz",
        snpsSummary="{dir}/snps_summary.txt.gz"
    shell:
        """
        gzip {input.snpsDepth}
        gzip {input.snpsFreq}
        gzip {input.snpsInfo}
        gzip {input.snpsSummary}
        """



# Combine species and genome information to create a taxonomy table
# for the species in the MIDAS database.
rule summarizeSpeciesTaxonomy:
    input:
        speciesInfo=config["speciesInfo"],
        genomeTaxonomy=config["genomeTaxonomy"]
    output:
        speciesTaxonomy=config["speciesTaxonomy"]
    script:
        "../scripts/runMIDAS/extractSpeciesTaxonomies.R"
