import pandas as pd 

rule insilicoseq:
    input:
        genome_file=config['genome_file'],
        abundance_file=lambda wildcards: df.loc[df['trial_parts'] == wildcards.trial_part,'abundance_file'],
    output:
        #workflow/out/insilicoSeqOutput/{trial_part}/insilico-reads.R1.fastq
        r1='workflow/out/insilicoSeqOutput20ones/{trial_part}/insilico-reads.fastq_R1.fastq',
        r2='workflow/out/insilicoSeqOutput20ones/{trial_part}/insilico-reads.fastq_R2.fastq'
    params:
        output_folder='workflow/out/insilicoSeqOutput20ones/{trial_part}/insilico-reads.fastq',
        total_reads = int(6666666),
    threads: config['maxCPUs']
    conda:
        "../../workflow/envs/insilico-no-builds-v3.yml"
    shell:
        """
        #mkdir {params.output_folder}
        iss generate --genomes {input.genome_file} --n_reads {params.total_reads} \
        --model novaseq --output {params.output_folder} --cpus {threads} \
        --abundance_file {input.abundance_file}
        """
