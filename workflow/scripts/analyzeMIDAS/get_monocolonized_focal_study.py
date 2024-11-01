import pandas as pd
import numpy as np
from glob import glob



def get_monocolonized(diversity_df, key_timepoint_df, subjects, species):
    div_initials = []
    subjects = []
    samples = []
    for i, subject in enumerate(subjects):
        good_samples = (key_timepoint_df.loc[key_timepoint_df['Subject'] == subject, 'Sample'].values)
        good_samples = [f'HouseholdTransmission-Stool-{sample}' for sample in good_samples]
        focal_samples = []
        focal_times = []
        for sample in good_samples:
            if sample in diversity_df['index'].values:
                focal_samples.append(sample)
                focal_times.append(int(sample.split('-')[-1]))
        print(focal_times)
        print(focal_samples)
        print(diversity_df['index'].values)
        if len(focal_times) == 0:
            continue
        min_sample = f'HouseholdTransmission-Stool-{subject}-{str(np.min(np.array(focal_times))).zfill(3)}'
        div_initial = diversity_df.loc[diversity_df['index'] == min_sample, '0']
        div_initials.append(div_initial)
        subjects.append(subject)
        samples.append(min_sample)
    monocolonized_df = pd.DataFrame(data = {'Div initial': div_initials, 'Subject': subjects, 'Sample': samples})
    monocolonized_df['Species'] = species
    return monocolonized_df


if __name__ == '__main__':

    fnames = glob('workflow/out/filtered_snps/*/diversity_df.csv')
    key_timepoint_df = pd.read_csv('tidysamplesKeyTimepoints.csv')

    monoclonized_initial_timepoint_dfs = []
    for fname in fnames:
        print(species)
        diversity_df = pd.read_csv(fname)
        subjects = list(diversity_df['Subject'].unique())
        print(subjects)
        species = list(diversity_df['Species'].unique())[0]
        monoclonized_initial_timepoint_df = get_monocolonized(diversity_df, key_timepoint_df, subjects, species)
        monoclonized_initial_timepoint_dfs.append(monoclonized_initial_timepoint_df)
    monoclonized_initial_timepoint_dfs = pd.concat(monoclonized_initial_timepoint_dfs )
    monoclonized_initial_timepoint_dfs.to_csv('initial_timepoint_diversity.csv')
    


    

    


        
    

