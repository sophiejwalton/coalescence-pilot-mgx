import glob 
import pandas as pd
import argparse



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Conduct Strain fishing for a given species across all cows')

    # add arguments
    parser.add_argument('minCoverage', action='store',
                    help='minCoverage for Species')


#    parser.add_argument('--main_study',
 #                   action='store_true')

    args = parser.parse_args()
    profiles = glob.glob('workflow/out/midasOutput/*/species/species_profile.txt')
    big_pd = []
    for fname in profiles:
        prof = pd.read_csv(fname, delimiter = '\t')
        sample = fname.split('/')[3]
        prof['Sample'] = sample
        prof['Subject'] = sample[:2]
        prof = prof.loc[prof['coverage'] > float(args.minCoverage), :]
        big_pd.append(prof)
    big_pd = pd.concat(big_pd)
    big_pd.to_csv("/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midasOutput/species/species_profile_all_abundant.csv")


