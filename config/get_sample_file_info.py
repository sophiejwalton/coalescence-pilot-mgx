import glob 
import pandas as pd



if __name__ == '__main__':
    sample_folders = glob.glob('/oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/*/')
    sample_names = [fname.split('/')[-2] for fname in sample_folders]
    read1 = [glob.glob(f'{folder}*R1*.gz')[0]for folder in sample_folders]
    read2 = [glob.glob(f'{folder}*R2*.gz')[0]for folder in sample_folders] 
    trim1= [f'/scratch/users/swalton/coalescence-pilot-mgx/trimmed/{sample}_L002-trimmed-pair1.fastq.gz' for sample in sample_names]
    trim2 = [f'/scratch/users/swalton/coalescence-pilot-mgx/trimmed/{sample}_L002-trimmed-pair2.fastq.gz' for sample in sample_names]
    df = pd.DataFrame(data = {'Sample': sample_names, 'Folder': sample_folders, 'read1': read1 , 'read2': read2, 'trim1': trim1, 'trim2': trim2})
    df.to_csv('sample_fnames.csv')
