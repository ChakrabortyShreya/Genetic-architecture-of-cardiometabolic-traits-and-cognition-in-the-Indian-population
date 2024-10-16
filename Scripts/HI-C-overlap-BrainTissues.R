dataset_name=c("GSE237345_ENCFF467KEZ_loops_GRCh38.bedpe","GSE237712_ENCFF730RTY_loops_GRCh38.bedpe","GSE238042_ENCFF216IAK_loops_GRCh38.bedpe","GSM5014501_HepG2_washU_text_Processed_Grch38.txt")
dataset_path="Datasets/"
gencode=read.table("gencode.v46.chr_patch_hapl_scaff.gene.processed.txt",header=TRUE)
str(gencode)
#GEO dictionary
GEO_accession=c("GSE237345_ENCFF467KEZ","GSE237712_ENCFF730RTY","GSE238042_ENCFF216IAK")
Tissue=c("Brain:Dorso-lateral_prefrontal_cortex","Brain:Caudate_nucleus","Brain:Posterior_cingulate_gyrus")
Gender=c("M","F","F")
Age=c(78,89,89)
data_dict=data.frame(GEO_accession,Tissue,Gender,Age)

#
library(tidyr)
library(dplyr)

#Function determines if credible set falls in bait or interacting region.
is_in_loop=function(a)
{
  #Subset by chromosome
  subset=hiC_D[hiC_D$BAIT==paste("chr",a[which(names(a)=="CHROM")],sep=""),]
  
  #Check credible set variant pos is in bait region or int region
  pos=as.numeric(a[which(names(a)=="POS")])
  pos_bait=subset[which(subset$BAIT_START<=pos & subset$BAIT_END>=pos),]
  pos_int=subset[which(subset$INT_START<=pos & subset$INT_END>=pos),]
  
  #If present in either region extract details
  if(dim(pos_bait)[1]!=0 | dim(pos_int)[1]!=0)
  {
    if(dim(pos_bait)[1]!=0 & dim(pos_int)[1]==0)
    {
      #To account for multiple bait/int regions
      n=dim(pos_bait)[1]
      s=cbind(do.call("rbind", replicate(n, a, simplify = FALSE)),pos_bait)
      return(s)
    }else if(dim(pos_bait)[1]==0 & dim(pos_int)[1]!=0)
    {
      #To account for multiple bait/int regions
      n=dim(pos_int)[1]
      s=cbind(do.call("rbind", replicate(n, a, simplify = FALSE)),pos_int)
      return(s)
    }else if (dim(pos_bait)[1]!=0 & dim(pos_int)[1]!=0)
    {
      #When bait/int overlaps
      n1=dim(pos_bait)[1]
      n2=dim(pos_int)[1]
      s1=cbind(do.call("rbind", replicate(n1, a, simplify = FALSE)),pos_bait)
      s2=cbind(do.call("rbind", replicate(n2, a, simplify = FALSE)),pos_int)
      s=unique(rbind(s1,s2))
      return(s)
    }
  }
}

#Fetch Gene TSS
fetch_gene_TSS=function(x)
{
  #If credible set is in bait region, fetch genes whose startsite-200bp falls in interacting region
  #print(x)
  if(x[which(names(x)=="Cred_Pos_BAIT")]=="YES")
  {
    start=as.numeric(x[which(names(x)=="INT_START")])
    end=as.numeric(x[which(names(x)=="INT_END")])
    chrom=paste("chr",x[which(names(x)=="CHROM")],sep="")
    get_gen=gencode[which(gencode$start_minus_200bp>=start &gencode$start_minus_200bp<=end & gencode$seqnames==chrom),]
    n_gen=dim(get_gen)[1]
    if(n_gen>0)
    {  
      gen_s=cbind(do.call("rbind", replicate(n_gen, x, simplify = FALSE)),get_gen)
      return(gen_s)
    }
  }else if (x[which(names(x)=="Cred_Pos_INT")]=="YES")
  {
    #If credible set is in int region, fetch genes whose startsite-200bp falls in bait region
    start=as.numeric(x[which(names(x)=="BAIT_START")])
    end=as.numeric(x[which(names(x)=="BAIT_END")])
    chrom=paste("chr",x[which(names(x)=="CHROM")],sep="")
    get_gen=gencode[which(gencode$start_minus_200bp>=start &gencode$start_minus_200bp<=end & gencode$seqnames==chrom),]
    n_gen=dim(get_gen)[1]
    if(n_gen>0)
    {
      gen_s=cbind(do.call("rbind", replicate(n_gen, x, simplify = FALSE)),get_gen)
      return(gen_s)
    }
  }
}

#Program starts here
phenotypes=c("MRT","Auditory_scaled","Visual_scaled","Dual_scaled",
             "Reading_scaled","Comprehension_scaled","Naming_Assoc_scaled",
             "Semantic","Phonetic","Vocab_scaled","ImmediateRecall_scaled",
             "DelayedRecall_scaled","NameRecog_scaled","NameFaceAssocFace_scaled",
             "Stroop_Level3_scaled","Span","Geo_Fig_scaled","TMT_B_A","HMSE","G_factor_scaled",
             "HDL","LDL","logTG","TC","TG_HDL","FBS","HbA1c","VisceralFat","VAI","MetSyn")
#phenotypes=c("HMSE")
length(phenotypes)

pheno_lg=c("HDL","LDL","logTG","TC","TG_HDL","FBS","HbA1c","VisceralFat","VAI","MetSyn")

cred_path="/Final_Credible_HiC_14052024_INDEX_HighConf/" //Path to clumped association summary statistics files

out_list=list()
for (j in 1:length(phenotypes))
{
  if (phenotypes[j] %in% pheno_lg)
  {
    cred_file_name=paste(phenotypes[j],".FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted_Final.txt",sep="")
  }else
  {
    cred_file_name=paste(phenotypes[j],"_LitCat.FamGrammarGamma.assoc_Processed_Indep_Hits_own.clumped_Formatted_Final.txt",sep="")
  }
  cred_file=read.table(paste(cred_path,cred_file_name,sep=""),sep="\t",header=TRUE)
  cred_file_v1= cred_file %>% separate(ID, into = c("CHROM", "POS","REF","ALT"), sep = ":", remove = FALSE)
  cred_file_v1$POS=as.numeric(cred_file_v1$POS)

  #Read HIC loop.bedpe
  out_cog=data.frame()
  for (i in 1:3)
  {
    hiC_D=read.table(paste(dataset_path,dataset_name[i],sep=""),sep="\t",header=TRUE)
    #if(length(which(hiC_D$BAIT_END>=hiC_D$INT_START))==0)
    #{
    #  print(paste("No Overlap between bait and interacting regions for GEO Accession:",i))
    #}
    new_data=apply(cred_file_v1,1,is_in_loop)
    if(!is.null(new_data))
    {
      new_data2 = new_data[-which(sapply(new_data, is.null))]
      #print(new_data2)
      new_data3=bind_rows(new_data2)
      #print(new_data3)
      new_data3$POS=as.numeric(new_data3$POS)
    
      new_data3$Cred_Pos_BAIT="NO"
      new_data3$Cred_Pos_BAIT[which((new_data3$POS>=new_data3$BAIT_START) & (new_data3$POS<=new_data3$BAIT_END))]="YES"
      #new_data3$Cred_Pos_BAIT[-which((new_data3$POS>=new_data3$BAIT_START) & (new_data3$POS<=new_data3$BAIT_END)]="NO"
      new_data3$Cred_Pos_INT=array()
      new_data3$Cred_Pos_INT[which(new_data3$POS>=new_data3$INT_START & new_data3$POS<=new_data3$INT_END)]="YES"
      new_data3$Cred_Pos_INT[-which(new_data3$POS>=new_data3$INT_START & new_data3$POS<=new_data3$INT_END)]="NO"
      new_data3$Cred_Pos_Overlap="NO"
      new_data3$Cred_Pos_Overlap[which((new_data3$Cred_Pos_BAIT=="YES") &(new_data3$Cred_Pos_INT=="YES"))]="YES"
      new_data3$GEO_accession=data_dict$GEO_accession[i]
      new_data3$Tissue=data_dict$Tissue[i]
      new_data3$Gender=data_dict$Gender[i]
      new_data3$Age=data_dict$Age[i]
      out_cog=rbind(out_cog,new_data3)
    }else
    {
      print(paste(phenotypes[j],":No overlap in Tissue",data_dict$Tissue[i]))
    }
  }
  out_cog_gen=apply(out_cog,1,fetch_gene_TSS)
  if(!is.null(out_cog_gen))
  {
    out_cog_gen_v2 = out_cog_gen[-which(sapply(out_cog_gen, is.null))]
    out_cog_gen_v3=bind_rows(out_cog_gen_v2)
    #Find out TSS of protein coding genes and the for those mapped to lead SNPs in our results
    out_cog_gen_v4=out_cog_gen_v3[((out_cog_gen_v3$gene_type=="protein_coding")&(out_cog_gen_v3$HighConf=="Yes")),]
    print("Hello")
    if(nrow(out_cog_gen_v4)>0)
    {
      out_cog_gen_v4$Phenotype=phenotypes[j]
      write.table(out_cog_gen_v4,paste("HIC-results/",phenotypes[j],"_hiC-Cog.txt",sep=""),quote=FALSE,row.names=FALSE)
      out_list[[j]]=out_cog_gen_v4
    }else
    {
      print(paste("Gencode v46, no protein coding genes mapped to loop for phenotype",phenotypes[j]))
    }
    
  }else
  {
    print(paste(phenotypes[j],":No protein coding TSS interacts"))
  }
}
Final_out=bind_rows(out_list)
write.table(Final_out,paste("HIC-results/All_results_hiC_Brain.txt",sep=""),quote=FALSE,row.names=FALSE)
