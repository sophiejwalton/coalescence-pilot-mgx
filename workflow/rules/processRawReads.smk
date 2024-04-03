# Use skewer to trim adapters from raw reads.
rule trimAdapters:
	input:
		r1=lambda wildcards: df.loc[str(wildcards.samplelane), 'read1'],
		r2=lambda wildcards: df.loc[str(wildcards.samplelane), 'read2'],
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
		trim1= lambda wildcards: df.loc[str(wildcards.samplelane), 'trim1'],
		trim2= lambda wildcards: df.loc[str(wildcards.samplelane), 'trim2'],
	conda:
		"../../workflow/envs/bowtie2-no-builds.yml"
	shell:
		"""
		bowtie2 --very-fast \
		  -x {params.humanref} \
		  -1 {params.trim1} -2 {params.trim2} \
		  -S {params.filterdir}/{wildcards.samplelane}.sam \
		  --un-conc-gz {params.filterdir}/{wildcards.samplelane}-filtered \
		  2> {output.bt2log}
		mv {params.filterdir}/{wildcards.samplelane}-filtered.1 {params.filterdir}/{wildcards.samplelane}-filtered.1.fastq.gz
		mv {params.filterdir}/{wildcards.samplelane}-filtered.2 {params.filterdir}/{wildcards.samplelane}-filtered.2.fastq.gz

		"""

#lambda wildcards: expand("workflow/out/midasOutput/species/abundantSpecies_{subject}.txt",
 #           subject=dict(tuple(df.groupby(['household'])))[wildcards.household]['subject'].tolist())

rule concatenate_fastq_files_across_lanes:
	input:
		filter1s=lambda wildcards: expand('workflow/out/filter/{wildcards.sample}-{lane}-filtered.1.fastq.gz',
			lane=[i.split('_')[-1][:4] for i in list(df.loc[df['Sample']==wildcards.sample,'filter1'].values)]), 
		filter2s=lambda wildcards: expand('workflow/out/filter/{wildcards.sample}-{lane}-filtered.2.fastq.gz',
			lane=[i.split('_')[-1][:4] for i in list(df.loc[df['Sample']==wildcards.sample,'filter2'].values)]), 
	output:
		concate1='workflow/out/concat/{sample}-filtered.1.fastq.gz',
		concate2='workflow/out/concat/{sample}-filtered.2.fastq.gz',
	threads: config['maxCPUs']
	params:
		concate_dir=config['concate_dir'],
		project=config['project'],
	shell:
		"""
		cat {input.filter1s} > {params.concate_dir}/{wildcards.sample}-filtered.1.fastq.gz
		cat {input.filter2s} > {params.concate_dir}/{wildcards.sample}-filtered.2.fastq.gz
		"""
