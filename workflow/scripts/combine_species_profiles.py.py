import glob 
import pandas as pd



if __name__ == '__main__':
    profiles = glob.glob('workflow/out/midasOutput/*/species/species_profile.txt')
    big_pd = []
    for fname in profiles:
        prof = pd.read_csv(fname)
        sample = fname.split('/')[3]
        prof['Sample'] = sample
        prof['Subject'] = sample[:2]
        big_pd.append(prof)
    big_pd = pd.concat(big_pd)
    big_pd.to_csv("/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midasOutput/species/species_profile_all.txt")

