# Use skewer to trim adapters from raw reads.
rule trimAdapters:
	input:
		r1=lambda wildcards: df.loc[str(wildcards.sample), 'read1'],
		r2=lambda wildcards: df.loc[str(wildcards.sample), 'read2'],
	params:
		outdir=config['trimdir'],
		project=config['project']
	output:
		trim1=join(config["trimdir"],"{sample}-trimmed-pair1.fastq.gz"),
		trim2=join(config["trimdir"],"{sample}-trimmed-pair2.fastq.gz"),
	threads: config['maxCPUs']
	conda:
		"../../workflow/envs/skewer-no-builds.yml"
	shell:
		"""
		skewer  -t {threads} -z -o {params.outdir}/{wildcards.sample} {input.r1} {input.r2}
		"""

# Use bowtie2 to filter out reads that map to the human genome.
# At this stage, also combine sequencing reads from all lanes.
rule filterOutHumanReads:
	input:
		trim1=lambda wildcards: df.loc[str(wildcards.samplelane), 'trim1'],
		trim2=lambda wildcards: df.loc[str(wildcards.samplelane), 'trim2']
	output:
		filtered1=join(config["filterdir"],"{samplelane}-filtered.1.fastq.gz"),
		filtered2=join(config["filterdir"],"{samplelane}-filtered.2.fastq.gz"),
		bt2log=join(config["filterdir"],"{samplelane}.bt2.log"),
		samfile=temp(join(config["filterdir"],"{samplelane}.sam"))
	threads: config['maxCPUs']
	params:
		humanref=config['humanGenomeRef'],
		filterdir=config['filterdir'],
		project=config['project'],
	#	trim1= lambda wildcards: df.loc[str(wildcards.samplelane), 'trim1'],
	#	trim2= lambda wildcards: df.loc[str(wildcards.samplelane), 'trim2'],
	conda:
		"../../workflow/envs/bowtie2-no-builds.yml"
	shell:
		"""
		bowtie2 --very-fast \
		  -x {params.humanref} \
		  -1 {input.trim1} -2 {input.trim2} \
		  -S {params.filterdir}/{wildcards.samplelane}.sam \
		  --un-conc-gz {params.filterdir}/{wildcards.samplelane}-filtered \
		  2> {output.bt2log}
		mv {params.filterdir}/{wildcards.samplelane}-filtered.1 {params.filterdir}/{wildcards.samplelane}-filtered.1.fastq.gz
		mv {params.filterdir}/{wildcards.samplelane}-filtered.2 {params.filterdir}/{wildcards.samplelane}-filtered.2.fastq.gz

		"""

#lambda wildcards: expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt",
 #           subject=dict(tuple(df.groupby(['household'])))[wildcards.household]['subject'].tolist())

