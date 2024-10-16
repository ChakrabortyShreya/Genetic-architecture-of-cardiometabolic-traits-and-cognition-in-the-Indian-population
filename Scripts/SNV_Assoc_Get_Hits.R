start=Sys.time()
print("Extracting association hits and QQplots......")
args<-commandArgs(TRUE)

#Enter path to Processed file
data_v1 = read.table(args[1], header=TRUE,sep=" ",fill=TRUE,row.names=NULL)
#head(data_v1)
str(data_v1)
print("Plotting QQ plots...")
#Genomic Inflation
lambda=qchisq(median(na.omit(data_v1$Pvalue)),1,lower.tail=FALSE)/qchisq(0.5,1)
print(paste("Genomic Inflation (lambda):",lambda))

#QQ_plot
data_obs=sort(na.omit(-log10(data_v1$Pvalue)))
N1=length(data_obs)

set.seed(18470)
Expected_data=sort(-log10(runif(N1,0,1)))
print("Expected_data")
length(Expected_data)

options(bitmapType='cairo')
tiff(paste(args[2],"/QQ_",args[3],".tiff",sep=""))
plot(Expected_data,data_obs,main=paste("QQ plot for",args[3]),ylab= "Observed quantiles of -log10(p)", xlab=" Empirical quantiles of -log10(p)",sub=paste(expression(lamdba),":",round(lambda,3)))
abline(coef=c(0,1))
dev.off()
print("Done")


#Get hits
#Get clumped file
data_v2= read.table(paste(args[1],"_Indep_Hits_own.clumped_Formatted_Final.txt_snpEff_Annotated.txt",sep=''),header=TRUE,fill=TRUE, sep='\t',row.names=NULL)
#str(data_v2)

GW_hits=data_v2[which(((data_v2$P < 5*(10^-8)) & (data_v2$Index.Clumped =="(INDEX)"))),]
write.table(GW_hits,paste(args[1],"Indep_GW_hits.txt",sep=""), quote=FALSE, row.names=FALSE)

#Suggestive hits at 5*10*-5
Sug_hits=data_v2[which(((data_v2$P < 5*(10^-5)) & (data_v2$Index.Clumped =="(INDEX)"))),]
write.table(Sug_hits,paste(args[1],"Indep_Sug_hits_5.txt",sep=""), quote=FALSE, row.names=FALSE)

#Suggestive hits at 5*10*-6
Sub_hits=data_v2[which(((data_v2$P < 5*(10^-6)) & (data_v2$Index.Clumped =="(INDEX)"))),]
write.table(Sub_hits,paste(args[1],"Indep_SubGenome_hits.txt",sep=""), quote=FALSE, row.names=FALSE)

#Suggestive hits at 5*10*-6 with SE<=6 and Info R2>0.7

Sub_hits_filter=Sub_hits[which((Sub_hits$INFO >0.7)&(Sub_hits$SE<=6)),]
write.table(Sub_hits_filter,paste(args[1],"Indep_SubGenome_SE6_INFO0.7_hits.txt",sep=""), quote=FALSE, row.names=FALSE)

#print("Hello")
#ManhattanPlot
#options(bitmapType='cairo')
library(qqman)
tiff(paste(args[2],"/ManhattanPlot_",args[[3]],".tiff",sep=""))
#print("tiff done")
high_snps=Sub_hits_filter$ID
manhattan(data_v1[which(is.na(data_v1$Pvalue)==FALSE),], chr="CHROM", snp="ID", bp="POS", p="Pvalue",suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8), highlight=high_snps, main=paste("Manhattan plot for ",args[[3]],sep=""))
#print("Here")
dev.off()
print("Done")
