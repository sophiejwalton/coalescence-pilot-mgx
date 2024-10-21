import glob 
import pandas as pd

def get_R2(x):
    x_list = x.split('_')
    x_list[-2] = 'R2'
    return '_'.join(x_list)


def get_trim1(x):
    string_sample = x.split('R1')[0][:-1]
    string_sample = string_sample.split('/')[-1]
    string_sample_well = string_sample.split('-')[0]
    if string_sample_well[0]=='A':
        string_sample_well = string_sample_well.replace('A','E')
    elif string_sample_well[0]=='B':
        string_sample_well = string_sample_well.replace('B','F')
    elif string_sample_well[0]=='C':
        string_sample_well = string_sample_well.replace('C','G')
    elif string_sample_well[0]=='D':
        string_sample_well = string_sample_well.replace('D','H')
    string_sample = string_sample_well + string_sample[2:]        
    
    return 'workflow/out/trimmed/' + string_sample + '-trimmed-pair1.fastq.gz'

def get_trim2(x):
    string_sample = x.split('trimmed-pair1')[0][:-1]
    string_sample = string_sample.split('/')[-1]
    
    return 'workflow/out/trimmed/' + string_sample + '-trimmed-pair2.fastq.gz'

def get_filter1(x):
    string_sample = x.split('trimmed-pair1')[0][:-1]
    string_sample = string_sample.split('/')[-1]
    
    return 'workflow/out/filter/' + string_sample + '-filtered.1.fastq.gz'

def get_filter2(x):
    string_sample = x.split('trimmed-pair1')[0][:-1]
    string_sample = string_sample.split('/')[-1]
    
    return 'workflow/out/filter/' + string_sample + '-filtered.2.fastq.gz'


#def get_sampleLane(x):
 #   filename = x.split('/')[-1]
   # return filename.split('_R1_001.fastq.gz')[0]
    

if __name__ == '__main__':

    R1 = glob.glob('/oak/stanford/groups/dpetrov/swalton/Sophie_Oct2024_Seq/*R1*.gz')
    df = pd.DataFrame(data = {'read1': R1})
    print(df)
    df['read2'] = df['read1'].transform(get_R2)
    df['trim1'] = df['read1'].transform(get_trim1)
    df['trim2'] = df['trim1'].transform(get_trim2)
    df['filter1'] = df['trim1'].transform(get_filter1)
    df['filter2'] = df['trim1'].transform(get_filter2)
   # df['SampleLane'] = df['read1'].transform(get_sampleLane)
    df['Sample'] = df['trim1'].transform(lambda x: x.split('-trimmed')[0][len('workflow/out/trimmed/'):])


    #df = pd.DataFrame(data = {'Sample': sample_names, 'read1': read1 , 'read2': read2, 'trim1': trim1, 'trim2': trim2})
    df.to_csv('sample_fnames_round2.csv')
