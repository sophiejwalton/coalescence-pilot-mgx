# Use skewer to trim adapters from raw reads.
rule trimAdapters:
	input:
		r1=lambda wildcards: df.loc[df['SampleLane'] == str(wildcards.samplelane), 'read1'],
		r2=lambda wildcards: df.loc[df['SampleLane'] == str(wildcards.samplelane), 'read2'],
	params:
		outdir=config['trimdir'],
		project=config['project']
	output:
		trim1=join(config["trimdir"],"{samplelane}-trimmed-pair1.fastq.gz"),
		trim2=join(config["trimdir"],"{samplelane}-trimmed-pair2.fastq.gz"),
	threads: config['maxCPUs']
	conda:
		"../../workflow/envs/skewer-no-builds.yml"
	shell:
		"""
		skewer  -t {threads} -z -o {params.outdir}/{wildcards.samplelane} {input.r1} {input.r2}
		"""

# Use bowtie2 to filter out reads that map to the human genome.
# At this stage, also combine sequencing reads from all lanes.
rule filterOutHumanReads:
	input:
		trim1=lambda wildcards: df.loc[wildcards.sample, 'trim1'],
		trim2=lambda wildcards: df.loc[wildcards.sample, 'trim2']
	output:
		filtered1=join(config["filterdir"],"{sample}-filtered.1.fastq.gz"),
		filtered2=join(config["filterdir"],"{sample}-filtered.2.fastq.gz"),
		bt2log=join(config["filterdir"],"{sample}.bt2.log"),
		samfile=temp(join(config["filterdir"],"{sample}.sam"))
	threads: config['maxCPUs']
	params:
		humanref=config['humanGenomeRef'],
		filterdir=config['filterdir'],
		project=config['project'],
		trim1= lambda wildcards: df.loc[wildcards.sample, 'trim1'],
		trim2= lambda wildcards: df.loc[wildcards.sample, 'trim2'],
	conda:
		"../../workflow/envs/bowtie2-no-builds.yml"
	shell:
		"""
		bowtie2 --very-fast \
		  -x {params.humanref} \
		  -1 {params.trim1} -2 {params.trim2} \
		  -S {params.filterdir}/{wildcards.sample}.sam \
		  --un-conc-gz {params.filterdir}/{wildcards.sample}-filtered \
		  2> {output.bt2log}
		mv {params.filterdir}/{wildcards.sample}-filtered.1 {params.filterdir}/{wildcards.sample}-filtered.1.fastq.gz
		mv {params.filterdir}/{wildcards.sample}-filtered.2 {params.filterdir}/{wildcards.sample}-filtered.2.fastq.gz
		"""
