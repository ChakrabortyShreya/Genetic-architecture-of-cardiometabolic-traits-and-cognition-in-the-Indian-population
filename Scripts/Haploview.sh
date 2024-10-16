#!/bin/bash

#Step1: Extract SV variants which are GW significant for lipid traits, which fall in genes GW significant in PTV associations (CETP,APOA5,APOE,TM6SF2,ZPR1)/ Similarly for AMIGO1
#Step2: remove indels
#Step3: Take variants only in genic regions of these 5 genes (NCBI coordinates)
#Step4: final_list in Final.txt


module load plink/2
module load plink/1.9.0
module load jdk-17.0.2

 
plink2 --extract CETP_ZPR1_APOA5_APOE_TM6SF2_extract_Final.txt --out Haplo_test_locus  --pfile ../../1.2.Imputation/ImputationResults/TLSA_SAS_MAC2/R2_gt_0.6_maf_0.0002/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink --make-bed

Rscript Dyslip_Fam_Update.R 

plink --bed Haplo_test_locus.bed --bim Haplo_test_locus.bim --fam Haplo_test_locus_update.fam --out Haplo_test_locus_pm --recode  --output-missing-phenotype 0


cut Haplo_test_locus_pm.map -f 2,4 > Haplo_test_locus_pm.info


java -jar ../Haploview.jar -nogui -pedfile Haplo_test_locus_pm.ped -info Haplo_test_locus_pm.info -assocCC -out Result -blockoutput ALL -log -dprime -svg

