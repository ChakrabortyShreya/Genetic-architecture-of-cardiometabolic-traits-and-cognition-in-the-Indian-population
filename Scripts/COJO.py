import os
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#Enter pheno
pheno=["HDL","LDL","TC","logTG","TG_HDL","VAI"]


#Enter File location
bfile_loc="/SANSCOG-2_108-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-maf0.0002_plink1.9"

if rank < len(pheno):
  #Extract independent significant hits
	os.system("plink --bfile "+bfile_loc+" --extract Conditional_Var/"+pheno[rank]+"_GW_reported.txt --indep 250 2 10 --out Conditional_Var/"+pheno[rank]+"_GW_Indep_VIF")
	os.system("gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --bfile "+bfile_loc+" --cojo-file COJO_FORMATTED_SUMSTATS/"+pheno[rank]+"_COJO_Form_Sumstats_25042024.ma --cojo-cond Conditional_Var/"+pheno[rank]+"_GW_Indep_VIF_25042024.prune.in --out Cond_Out/VIF_Filtered/"+pheno[rank]+"_GW_cond_VIF_25042024 --threads 10")
  
