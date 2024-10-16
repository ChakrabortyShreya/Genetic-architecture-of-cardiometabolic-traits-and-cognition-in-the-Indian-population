library(enrichR)
websiteLive <- getOption("enrichR.live")
#if (websiteLive) {
#   listEnrichrSites()
#    setEnrichrSite("Enrichr") # Human genes   
#}

#if (websiteLive) dbs <- listEnrichrDbs()
#if (websiteLive) dbs

args = commandArgs(trailingOnly=TRUE)
input_genes=unique(read.table(args[[1]],header=FALSE)$V1)
phenotype=args[[2]]
set_number=args[[3]]

library(plyr)

#Specify database
#to_use_dbs <- c("KEGG_2021_Human","GeDiPNet_2023","WikiPathway_2023_Human","huMap","SYNGO_2022","OMIM_Expanded")

to_use_dbs <- c("KEGG_2021_Human","Reactome_2022","SynGO_2024")

#print("Hello")
plot_top20_pathways=function(db,pheno,set)
{
	library(ggplot2)
	library(stringr)
	db$Gene_count=count.fields(textConnection(db$Genes), sep = ";")
	str(db)
	sorted_data=db[order(db$Odds.Ratio, decreasing = TRUE), ] 
	#if (nrow(sorted_data)> 20)
	#{
	#	sorted_data=sorted_data[c(1:20),]
	#}
	#term_order=sorted_data$Term[order(sorted_data$Adjusted.P.value,sorted_data$.id)]
	#sorted_data$Term=factor(sorted_data$Term, levels = unique(term_order))
	
  sorted_data$.id <- factor(sorted_data$.id , levels = to_use_dbs, ordered = TRUE)
	sorted_data_subset <-sorted_data %>% group_by(.id) %>% slice_head(n = 10) 
	write.table(sorted_data_subset,paste("Enriched/",pheno,"_enrichr.txt",sep=""),row.names=FALSE)
	
	ggplot(sorted_data_subset, aes(x = Gene_count, y = Term, fill = Adjusted.P.value)) +
  		geom_bar(stat = "identity") + 
  		scale_fill_viridis_c() +
		theme_bw() +
		scale_colour_brewer(palette = "Dark2") +
		theme(panel.grid.major.y = element_blank()) + 
		labs(x="Gene Count", y="Enriched Term")+
		facet_grid(.id ~ ., scales = "free_y", space="free_y") +
		#scale_x_continuous(breaks= c(1:(max(sorted_data$Gene_count)+3)))+
		scale_x_continuous(breaks= c(1:(max(sorted_data$Gene_count))))+
		scale_y_discrete(label = function(x) str_wrap(str_trunc(x, 40), width = 30))+
		ggtitle(paste("Top 15 enriched pathways for ",pheno," (",set,")",sep=""))
	ggsave(paste("Plots/",set,"/",pheno,"_",set,"_enrichr.tiff",sep=""), units="in", dpi=300, height=7, width =10,compression = 'lzw')
	#dev.off()
	
}

#====================Main Enricher call ===========================
library(dplyr)
#Function to run Enrichr across all databases
run_enrichr=function(gene_list,pheno,set)
{
	print("Hello!")
	if (websiteLive) {
		enriched <- enrichr(gene_list, to_use_dbs)
    		enriched_all <- ldply(enriched, data.frame)
		print(head(enriched_all))
		sig_enriched=enriched_all[enriched_all$Adjusted.P.value < 0.05,]
		write.table(sig_enriched,"test.txt",row.names=FALSE, quote=FALSE,sep="|")
		print(paste("Enrichment results obtained for ",pheno," - ",set,sep=""))
		print("Plotting Results...")
		plot_top20_pathways(sig_enriched,pheno,set)	#Call plot function
		print("Done.")
	}
}

run_enrichr(input_genes,phenotype,set_number)
#===================================================================
