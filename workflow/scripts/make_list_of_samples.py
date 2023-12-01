import glob 
import pandas as pd
import argparse



if __name__ == '__main__':
    good_samples_folders = glob.glob('/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midas2_output/*/species/species_profile.tsv')
    samples = [sample_folder.split('/')[10] for sample_folder in good_samples_folders]
    df = pd.DataFrame(data = {'sample_name': samples}).set_index('sample_name')
    df['midas_outdir'] = '/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/midas2_output/'
    df.to_csv('/oak/stanford/groups/dpetrov/swalton/CZpilot-seq-invitro/workflow/out/list_of_samples.tsv', sep="\t")
