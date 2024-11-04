import glob
import pandas as pd
import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Get changing SNPS within populations')
  parser.add_argument('save_Dir', action = 'store')
  parser.add_argument('fname_find', action = 'store')
  parser.add_argument('fname_save', action = 'store')
  args = parser.parse_args()
  fnames = glob.glob(f'{args.save_Dir}/*/{args.fname_find}')
  big_df = []
  for fname in fnames:
    small_df = pd.read_csv(fname)
    if len(small_df) > 0:
      big_df.append(pd.read_csv(fname))
      
  big_df = pd.concat(big_df)
  big_df.to_csv(args.fname_save)
