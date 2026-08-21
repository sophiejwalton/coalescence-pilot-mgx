import pandas as pd
import gzip 
import numpy as np 

import gzip
import time
# google answer
def num_reads(filename):
    count = 0
    with gzip.open(filename, 'rt') as f:
        for i, _ in enumerate(f):
            pass
    return (i + 1) // 4


if __name__ == '__main__':
    df = pd.read_csv('TableS1_only1.csv')
    start = time.time()
    df['fname2'] = df['fname'].transform(lambda x: x.split('.1.fastq.gz')[0]+ '.2.fastq.gz')
    df['accession'] = 'TBA'
    df['num_reads1'] = np.nan
   # df['num_reads2'] = np.nan
    df = df.set_index('fname')
    for i, fname in enumerate(df.index.values[:10]):
        print(i, df.loc[fname,'fname2'])
        df.loc[fname,'num_reads1'] = num_reads(fname)
       # df.loc[fname,'num_reads2'] = num_reads(df.loc[fname,'fname2'])
       # print(df.loc[fname,:])
  #  df['total_reads'] = df['num_reads1'] + df['num_reads2'] 
    df.to_csv('TableS1_with_num_reads_only1.csv')
    print(time.time()-start)
  
