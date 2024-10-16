
args=commandArgs(trailingOnly=TRUE)


#Instruments=read.table(args[1],header=TRUE)
#Outcome=read.table(args[2],header=TRUE)

Pheno=args[2]

library(TwoSampleMR)
library(ggplot2)


methods=mr_method_list()
#methods$obj

print("Reading Instruments....")

# These columns should be present in the exposure GWAS file
exposure_dat=read_exposure_data(args[1],clump = FALSE,
  sep = "\t",
  phenotype_col = "Pheno",
  snp_col = "ID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  samplesize_col = "N",
  gene_col = "Gene Name",
  id_col = "ID",
  log_pval = FALSE,
  chr_col = "CHROM",
  pos_col = "POS"
)

#Reads outcomes and stores in a compatible format
print("Reading Outcome...")
outcome=read.table(args[2],header=TRUE)
colnames(outcome)=c("chr.outcome","pos.outcome","other_allele.outcome", "effect_allele.outcome","samplesize.outcome","eaf.outcome","beta.outcome","betavar.outcome", "pval.outcome","SNP","se.outcome")


#Filters outcome sumstats for Instruments
outcome$SNP=tolower(outcome$SNP)
print("Extracting outcome summary statistics for instruments....")
outcome_dat=outcome[which(outcome$SNP %in% exposure_dat$SNP),]
outcome_dat$outcome=rep(args[3],nrow(outcome_dat))

#Required columns added
print("Assigning outcome ID...")
library(dplyr)
outcome_dat$id.outcome=rep(paste(exposure_dat$id.exposure[1],"_o",sep=""),nrow(outcome_dat))

print(outcome_dat)

#Harmonize the data
har_dat=harmonise_data(exposure_dat, outcome_dat, action = 1)
print("Harmonized data.")
print(har_dat)


out_name=paste(exposure_dat$exposure[1],".",outcome_dat$outcome[1],sep="")
#print(out_name)

#MR analysis
print("All MR Results.....")
#res_all=mr(har_dat, method_list=methods$obj)
#res_all=mr_wrapper(har_dat)
#res_all
#str(res_all)
#print("Hello")
id_name=paste(exposure_dat$id.exposure[1],".",outcome_dat$id.outcome[1],sep="")
print(id_name)

#get_res=res_all[[id_name]]
#MR_est=get_res[["estimates"]]
#MR_est$Exposure=rep(exposure_dat$exposure[1],nrow(MR_est))
#MR_est$Outcome=rep(outcome_dat$outcome[1],nrow(MR_est))

#print(get_res)

library(ggplot2)
res=mr(har_dat,method_list=c("mr_ivw","mr_egger_regression","mr_ivw_fe","mr_ivw_mre"))
res
get_res=res[[id_name]]
MR_est=get_res[["estimates"]]
#MR_est
#MR_est$Exposure=rep(exposure_dat$exposure[1],nrow(MR_est))
#MR_est$Outcome=rep(outcome_dat$outcome[1],nrow(MR_est))


p1 <- mr_scatter_plot(res, har_dat)
ggsave(p1[[1]], file = paste("Plots_25042024/",out_name,"_Effects.png",sep=""), width = 7, height = 7,dpi=300,device="png")

res_single <- mr_singlesnp(har_dat,single_method = "mr_wald_ratio",all_method = c("mr_ivw", "mr_egger_regression","mr_ivw_fe","mr_ivw_mre"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = paste("Plots_25042024/Forest_plot_",out_name,"_Effects.png",sep=""), width = 7, height = 7,dpi=300,device="png")


p3 <- mr_funnel_plot(res_single)
ggsave(p3[[1]], file = paste("Plots_25042024/Funnel_plot_",out_name,"_Effects.png",sep=""), width = 7, height = 7,dpi=300,device="png")

library(gridExtra)
combined_plot <- arrangeGrob(p1[[1]], p3[[1]], nrow = 1)
ggsave(combined_plot, file=paste("Plots_25042024/MR_Plots_",out_name,".png",sep=""),width = 7, height = 7,dpi=300,device="png")


#Heterogeneity Test
print("Conducting heterogeneity Tests...")
het=mr_heterogeneity(har_dat,method_list =c("mr_two_sample_ml","mr_egger_regression","mr_ivw","mr_uwr","mr_ivw_fe","mr_ivw_mre"))
het

#test horizontal pleiotropy
print("Horizontal Pleiotropy Results...")
hor_plei=mr_pleiotropy_test(har_dat)
str(hor_plei)
hor_plei


#Single test Forest plot
#res_single <- mr_singlesnp(har_dat)
#p2 <- mr_forest_plot(res_single)
#ggsave(p2[[1]], file = paste("Plots/",outname,"_SingleVsAll_FP.png"), width = 7, height = 7)

out=directionality_test(har_dat)

library(openxlsx)
wb = createWorkbook()
addWorksheet(wb, "MR_All")
addWorksheet(wb, "Heterogeneity")
addWorksheet(wb, "Horizontal_pleiotropy")
addWorksheet(wb, "Directionality")
writeData(wb, sheet = "MR_All", x = res)
writeData(wb, sheet = "Heterogeneity", x = het)
writeData(wb, sheet = "Horizontal_pleiotropy", x = hor_plei)
writeData(wb, sheet = "Directionality", x = out)
saveWorkbook(wb, paste("Results_25042024/MR_results_",out_name,".xlsx",sep=""),overwrite=TRUE)
