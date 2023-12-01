# Analyze the SNP calls for each species.
# Process the SNP calls, plot the SFS,
# and analyze strain similarity between subjects
# using both fixed differences and strain fishing methods.
rule filterVariantsCalculateStrainSimilarity:
    input:
        snpsDepth="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_depth.txt.gz",
        snpsFreq="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_freq.txt.gz",
        snpsInfo="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_info.txt.gz",
        snpsSummary="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/{species}/snps_summary.txt.gz"
    output:
        "workflow/report/plotMIDAS/calculateStrainSimilarity/{species}/done.txt"
    params:
        plotDefaults=config["plotDefaults"],
        indir="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",
        outdir="workflow/out/midasOutput/snps/HouseholdTransmission-Stool/",
        plotdir="workflow/report/plotMIDAS/calculateStrainSimilarity/",
        species="{species}",
        runAll="TRUE"
    conda:
        "../../workflow/envs/Renv-minimal.yml"
    script:
        "../scripts/plotMIDAS/filterVariantsCalculateStrainSimilarity.R"

# Use the species relative abundances calculated by MIDAS
# to plot changes in community composition over time.
# rule plotSpeciesAbundances:
#    input:
#        speciesAbundances="workflow/out/midasOutput/species/species_profile_all.txt",
#        speciesTaxonomy=config["speciesTaxonomy"]
#    output:
#        plotpdf=join(config["plotMIDASoutdir"],"sampleSpeciesAbundances.pdf"),
#        plotpng=join(config["plotMIDASoutdir"],"sampleSpeciesAbundances.png"),
#        plotLogpdf=join(config["plotMIDASoutdir"],"sampleSpeciesAbundances-log.pdf"),
#        plotLogpng=join(config["plotMIDASoutdir"],"sampleSpeciesAbundances-log.png"),
#        plotRelativeAbundanceByRankpng=join(config["plotMIDASoutdir"],"sampleRelativeAbundanceByRank.png")
#    params:
#        plotDefaults=config["plotDefaults"]
#    script:
#        "../scripts/plotMIDAS/plotSpeciesAbundances.R"
