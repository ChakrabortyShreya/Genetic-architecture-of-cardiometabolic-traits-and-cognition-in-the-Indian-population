import os
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#Phenotypes to test
z=["MRT","Auditory_scaled","Visual_scaled","Dual_scaled","Reading_scaled","Comprehension_scaled","Naming_Assoc_scaled","Semantic","Phonetic","Vocab_scaled","ImmediateRecall_scaled","DelayedRecall_scaled","NameRecog_scaled","NameFaceAssocName_scaled","NameFaceAssocFace_scaled","Stroop_Level3_scaled","Span","Geo_Fig_scaled","Matrix_raw"]


if rank < len(z):
	
	os.system("mkdir "+z[rank]+"\n")
	with open(z[rank]+"/"+z[rank]+".log","w") as f:
		f.write("Running rvtests...\n")
		os.system("time rvtest --inVcf SANSCOG-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002.vcf.gz  --out "+z[rank]+"/2_108_"+z[rank]+"_LitCat_corrected 
              --covar ../../1.4.Phenotypes/MAC2_R2_MAF0.0002/COGNITO_EduDum_pheno_MAC2_R2_MAF0.0002_pca.txt  --covar-name GEN_SEX,Age_FINAL_COGNITO,Age_FINAL_COGNITO_sq,LitMed,LitHigh,U1,U2,U3,U4,U5
              --pheno ../../1.4.Phenotypes/MAC2_R2_MAF0.0002/COGNITO_EduDum_pheno_MAC2_R2_MAF0.0002_pca.txt --inverseNormal --pheno-name "+z[rank]+" 
              --single famGrammarGamma --kinship ../../1.3.Relatedness/TLSA_Imp_MAC2_R2_0.6_maf_0.0002/SANSCOG-2_108_Imp_MAC2_R2_0.6_MAF_0.0002_pruned_bn.kinship 
              --numThread 15 --useResidualAsPhenotype")
		f.write("Done.\n")
		f.write("Processing GWAS sumstats..(Calc Z SE)...\n")
		os.system("time Rscript ../1.Clump_Pre_process.R  "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc")
		f.write("Done.\n")
		f.write("Performing LD clump...\n")
		os.system("time plink --bfile SANSCOG-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9 --clump "+z[rank]+"/2_108_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed 
              --clump-best --clump-field Pvalue --clump-kb 250 --clump-p1 5e-6 --clump-p2 5e-5 --clump-r2 0.7 --clump-snp-field ID --clump-verbose 
              --out "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own  --threads 30")
		f.write("Done.\n")
		f.write("Processing .clumped file....\n")
		os.system("time python ../2.Clump_post_process.py "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped > "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted.txt")
		f.write("Done.\n")
		f.write("Getting Independent hits and clumps....\n")
		os.system("time python ../3.Get_Indep.py "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted.txt "+z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed "+ z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted_Final.txt")
		f.write("Done.\n")
		f.write("GWAS summary statistics is ready!\n")

		f.write("Fetching SNPEFF annotations on Independent clumps....\n")
		os.system("time python ../4.SNP_EFF_Annot.py "+ z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted_Final.txt")
		f.write("Done.\n")
		f.write("Writing Tables and Plotting QQ plots...\n")
		os.system("time Rscript ../5.Get_hits.R "+ z[rank]+"/2_108_NoAPOE_"+z[rank]+"_LitCat_corrected.FamGrammarGamma.assoc_Processed "+z[rank]+"/ "+z[rank])
		f.write("Done.")
		f.close()
