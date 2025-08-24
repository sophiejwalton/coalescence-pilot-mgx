#!/bin/bash

#SBATCH --job-name=testgt-pro
#SBATCH --time=1440
#SBATCH -p normal
#SBATCH -c 4

ml system rclone

rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/AG-Fecal_S551_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/AG-Fecal_S551_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/AM-Fecal_S552_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/AM-Fecal_S552_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BA-Fecal_S553_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BA-Fecal_S553_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BB-Fecal_S554_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BB-Fecal_S554_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BC-Fecal_S555_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BC-Fecal_S555_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BF-Fecal_S556_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BF-Fecal_S556_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BG-Fecal_S558_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BG-Fecal_S558_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BK-Fecal_S557_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/data/coalescence-pilot/coalescence_data/BK-Fecal_S557_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AA_S365_L003_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AK_S369_L003_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AA_S365_L003_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AK_S369_L003_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AA_S365_L004_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AK_S369_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AA_S365_L004_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AK_S369_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AC_S366_L003_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AR_S370_L003_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AC_S366_L003_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AR_S370_L003_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AC_S366_L004_R1_001.fastq.gz gdrive:sophie_data 
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AR_S370_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AC_S366_L004_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-AR_S370_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AE_S367_L003_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BD_S371_L003_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AE_S367_L003_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BD_S371_L003_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AE_S367_L004_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BD_S371_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AE_S367_L004_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BD_S371_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AF_S368_L003_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BE_S372_L003_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AF_S368_L003_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BE_S372_L003_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AF_S368_L004_R1_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BE_S372_L004_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e003-inoculumn-fecal-AF_S368_L004_R2_001.fastq.gz gdrive:sophie_data  
rclone copy /oak/stanford/groups/dpetrov/Feb2024_biohub_SJW_JCCV_KAS/Sophie_Walton/e004-inoculumn-fecal-BE_S372_L004_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/Feb2025_Seq/Assembly-G12-fecal-AA-fecal-0_S703_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/Feb2025_Seq/Coalescence-H12-AA-fecal-0_S715_R1_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/Feb2025_Seq/Assembly-G12-fecal-AA-fecal-0_S703_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/Feb2025_Seq/Coalescence-H12-AA-fecal-0_S715_R2_001.fastq.gz gdrive:sophie_data
rclone copy /oak/stanford/groups/dpetrov/swalton/CZSeq

