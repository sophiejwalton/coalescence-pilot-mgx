import glob 
import pandas as pd
import argparse
import os 


if __name__ == '__main__':
    fnames = glob('*/*snps_*.tsv.lz4')
    for fname in fnames:
        command = f'lz4 {fname}'
        os.system(command)
    


   
