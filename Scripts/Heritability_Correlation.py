import os

bfile_loc="/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9"

#Segment-based LD score
#os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile "+bfile_loc+" --ld-score-region  200 --ld-wind 10 --ld-rsq-cutoff 0.7 --out LD_seg_score_200_10_0.7 --threads 60")
#os.system("Rscript  Stratify_Seg_LD_score.R LD_seg_score_200_10_0.7.score.ld")

#Creating GRM matrix
for i in ["1","2","3","4"]:
	os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile "+bfile_loc+" --extract GRM/snp_group"+i+".txt  --make-grm --out GRM/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9_group"+i+"_GRM --threads 60")


#Heritability
pheno1=["HDL","LDL","TG_HDL","logTG","TG","TC","FBS","HbA1c"] #Similarly cognitive phenotypes
for i in range(len(pheno1)):
	os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --reml --mgrm GRM/multi_GRMs.txt --pheno INV_res_pheno/INV_cov_adj_"+pheno1[i]+"_pheno_no_fid.txt  --out GREML_Her/GREML_Her_"+pheno1[i]+"_25042024 --threads 60")
 
#Genetic correlation

os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --mgrm GRM/multi_GRMs.txt --make-grm  --out GRM/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9_GRM --threads 60")

for i in range(0,len(pheno)):
	for j in range(0,len(pheno)):
		if i<j:
			os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --reml-bivar --grm GRM/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9_GRM --pheno INV_res_pheno_pair/INV_cov_adj_"+pheno[i]+"_"+pheno[j]+"_no_fid.txt  --out GREML_Cor/GREML_Cor_"+pheno[i]+"_"+pheno[j]+" --threads 60")

