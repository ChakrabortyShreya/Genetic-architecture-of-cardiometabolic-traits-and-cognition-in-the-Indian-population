import os
import sys
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

set_entered=sys.argv[1]

#Enter folder paths here:
base_dir=""  
imp_file=f'{base_dir}/SANSCOG-imputed-filtered-r0.6-sorted-mac2-annot-atgc-rem-maf0.0002.vcf.gz'

sanscog_bed_dir=f'{base_dir}/2.0.Gene_Based_Assoc/REGENIE_bed_files/'
sanscog_bed_fullset=f'SANSCOG-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9'
sanscog_bed_subset=f'SANSCOG-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9_MAF0.01'

phen_dir=f'{base_dir}/Pheno_Covar/'
out_dir=f'{base_dir}/Out/'
log_dir=f'{base_dir}/Log/'


#Enter Testing categories
#sets=["coding","PTV"]


#Create output folders
#os.system(f'mkdir {out_dir}/Coding {out_dir}/PTV')

#Enter Phenotypes
phenotypes=["MRT", "Auditory_scaled", "Visual_scaled", "Dual_scaled", "Reading_scaled", "Comprehension_scaled", "Naming_Assoc_scaled", "Semantic", "Phonetic", "Vocab_scaled", "ImmediateRecall_scaled", "DelayedRecall_scaled", "NameRecog_scaled", "NameFaceAssocName_scaled", "NameFaceAssocFace_scaled", "Stroop_Level3_scaled", "Span", "Geo_Fig_scaled", "Matrix_raw","G_factor_scaled","BMI", "WHR", "WHeightR", "BodyFat", "VisceralFat", "VAI", "WHRAdjBMI", "FBS", "HbA1c", "TC", "TG", "HDL", "LDL","logTG", "TG_HDL", "G_factor_scaled","TMT_B_A","Dyslipidemia","HMSE"]


#Function to run regenie with phenotype and set supplied as arguments
def analysis(pheno,set_number):
#	print(pheno+"..."+set_number)
	os.system(f'export LC_ALL="en_US.UTF-8"')
	anno_file=f'{base_dir}/Final_SNPEFF_Annot/sanscog_annotations_'+set_number+'.txt.gz'
	set_file=f'{base_dir}/Final_SNPEFF_Annot/sanscog_annotations_'+set_number+'_set_list.txt.gz'
	mask_file=f'{base_dir}/Final_SNPEFF_Annot/snpEff_mask_definitions_consequence_'+set_number+'.txt.gz'
	variant_ids=f'{base_dir}/Final_SNPEFF_Annot/variantIDs_'+set_number+'.txt'
	
	if pheno!='Dyslipidemia':
		#REGENIE Step1
		os.system(f"time regenie --step 1 --bed "+sanscog_bed_dir+sanscog_bed_subset+" --phenoFile "+phen_dir+"cov_adj_"+pheno+"_pheno_no_fid.txt  --phenoCol cov_adj_"+pheno+" --apply-rint --force-step1 --bsize 1000 --qt --lowmem --extract "+variant_ids+" --out "+out_dir+set_number+"/"+pheno+"_step1_"+set_number+"_output")
		#REGENIE Step2
		os.system(f"time regenie --step 2 --bed "+sanscog_bed_dir+sanscog_bed_fullset+" --phenoFile "+phen_dir+"cov_adj_"+pheno+"_pheno_no_fid.txt --phenoCol cov_adj_"+pheno+" --apply-rint --anno-file "+anno_file+"  --set-list "+set_file+"  --mask-def "+mask_file+" --extract "+variant_ids+"  --bsize 1000  --qt  --approx --pThresh 0.01   --pred "+out_dir+set_number+"/"+pheno+"_step1_"+set_number+"_output_pred.list  --aaf-bins 0.001,0.01,0.05 --vc-tests acato,skato-acat  --out "+out_dir+set_number+"/"+pheno+"_step2_"+set_number+"_output")	

#	else:
		#REGENIE Step1-Binary Phenotypes		
		#os.system(f"regenie   --step 1  --bed {sanscog_bed_dir}/{data_bed_subset_file_stem}  --phenoFile {pheno_file} --covarFile {cov_file} --force-step1 --bsize 1000 --bt --lowmem --extract {variant_ids}  --out {phen_dir}/{entry}_step1_set1_output &> {phen_dir}/{entry}_step2_set2_output.log")
		
		#REGENIE Step2-Binary Phenotypes   
#		os.system(f"time regenie  --step 2  --bed {sanscog_bed_dir}/sanscog_2_108_bed_fullset --covarFile {cov_file} --phenoFile {pheno_file} --anno-file {anno_file}  --set-list {set_file}  --mask-def {mask_file} --extract {variant_ids}  --bsize 1000  --bt  --approx --pThresh 0.01   --pred {phen_dir}/{entry}_step1_set2_output_pred.list  --aaf-bins 0.001,0.01,0.05 --vc-tests acato,skato-acat  --out {phen_dir}/{entry}_step2_set2_output &> {phen_dir}/{entry}_step2_set2_output.log")



#Preprocess for parallelization across nodes and cores
#pheno_file_parallel= [ele for ele in phenotypes for i in range(3)]
#set_number=sets*len(phenotypes)


#def test(pheno,set_number):
#	line=pheno+"..."+set_number
#	return line

#analysis(pheno_file_parallel[0],set_number[2])

#for rank in range(len(pheno_file_parallel)):
for rank in range(len(phenotypes)):
	
	log_path = os.path.join(log_dir, f"{phenotypes[rank]}_{set_entered}_Regenie_RunTime.log")	
	log_file = open(log_path, 'w')
	log_file.write("Gene-based analysis started for "+phenotypes[rank])
	log_file.write("Set="+set_entered)
	#line2=test(pheno_file_parallel[rank],set_number[rank])
	#log_file.write(line2)
	start_time = time.time()
	#FUNCTION CALL
	analysis(phenotypes[rank],set_entered)
	end_time = time.time()
	elapsed_time = round((float((end_time - start_time)/60)),3)
	line=f"Regenie Analysis  \n \ttook {elapsed_time} minutes\n"
	log_file.write(line)
	log_file.write("Analysis completed. Have a nice day!")
