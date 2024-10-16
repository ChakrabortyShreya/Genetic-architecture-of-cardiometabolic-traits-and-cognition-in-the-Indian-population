library(coloc)

# Specify the input directory containing the dataset1 files
input_directory <- "Input_datasets/"

# Specify the pattern for dataset1 and dataset2 file names
dataset1_pattern <- "common_(.*?)(_dataset1\\.csv)"  #GWAS SUMMARY STATS
dataset2_pattern <- "common_\\1_dataset2.csv"        #EQTL NES, SD

# Specify the output directory
output_directory <- "COLOC_BETA_MAF_N_eQTL_BETAVAR_SDY/"

# Get a list of dataset1 files matching the pattern
dataset1_files <- list.files(path = input_directory, pattern = dataset1_pattern, full.names = TRUE)
#
# Create the output directory if it doesn't exist
dir.create(output_directory, showWarnings = FALSE)

# Loop through each dataset1 file
for (dataset1_file in dataset1_files) {
  # Extract dataset prefix from the file name
  dataset_prefix <- sub(dataset1_pattern, "\\1", basename(dataset1_file))

  # Construct the corresponding dataset2 file name
  dataset2_file <- file.path(input_directory, sub(dataset1_pattern, dataset2_pattern, basename(dataset1_file)))

  # Check if the corresponding dataset2 file exists
  if (file.exists(dataset2_file)) {
    # Load dataset1
    dataset1 <- read.csv(dataset1_file, header = TRUE)

    # Load dataset2
    dataset2 <- read.csv(dataset2_file, header = TRUE)

    # Get the list of unique common SNPs
    common_snps <- intersect(unique(dataset1$rsid), unique(dataset2$rsid))

    # Create an empty list to store results for the current combination
    results_list <- list()

    # Loop through each common SNP
    for (rsid1 in common_snps) {
      # Subset dataset1 for the current rsid
      subset_dataset1 <- subset(dataset1, rsid == rsid1)

      # Subset dataset2 for the current rsid
      subset_dataset2 <- subset(dataset2, rsid == rsid1)
#
      # Check if common SNPs are present in both datasets
      if (nrow(subset_dataset1) > 0 && nrow(subset_dataset2) > 0) {
        # Run coloc.abf for each row with unique SNP
        for (row_index in 1:nrow(subset_dataset2)) {
          # Extract information for the current row
          Y1 <- subset_dataset1$Beta
          Y2 <- subset_dataset2$nes[row_index]
          N1 <- subset_dataset1$N
          N2 <- subset_dataset2$n[row_index]
          p1 <- as.numeric(subset_dataset1$pval)
          p2 <- as.numeric(subset_dataset2$pval[row_index])
          MAF1 <- subset_dataset1$MAF
          MAF2 <- subset_dataset2$maf[row_index]
          snp1 <- subset_dataset1$rsid
          snp2 <- subset_dataset2$rsid[row_index]
          #se1 <- subset_dataset1$sdY
          #se2 <- subset_dataset2$se[row_index]
	  V1 <- subset_dataset1$Varbeta
	  V2 <- subset_dataset2$betavar[row_index]
          geneSymbol <- subset_dataset2$geneSymbol[row_index]
          tissue <- subset_dataset2$tissue[row_index]
          Index.Clumped <- subset_dataset2$Index.Clumped[row_index]
          CHROM_y <- subset_dataset1$CHROM_y
          POS_y <- subset_dataset1$POS_y
          REF_y <- subset_dataset1$REF_y
          ALT_y <- subset_dataset1$ALT_y
          Index.Clumped <- subset_dataset1$Index.Clumped
          INDEX_ID <- subset_dataset1$Existing_Variation
          CHROM_x <- subset_dataset2$CHROM_x[row_index]
          POS_x <- subset_dataset2$POS_x[row_index]
          REF_x <- subset_dataset2$REF_x[row_index]
          ALT_x <- subset_dataset2$ALT_x[row_index]
          se <- subset_dataset2$se[row_index]
          nes <- subset_dataset2$nes[row_index]

          # Check for missing values
          
	   #if (any(is.na(Y1)) || any(is.na(Y2)) || any(is.na(N1)) || any(is.na(N2)) || any(is.na(p1)) || any(is.na(p2)) ||
           #	any(is.na(MAF1)) || any(is.na(MAF2)) || any(is.na(se1)) || any(is.na(se2))) {
	    if (any(is.na(Y1)) || any(is.na(Y2)) || any(is.na(N1)) || any(is.na(N2)) || any(is.na(p1)) || any(is.na(p2)) ||
              any(is.na(MAF1)) || any(is.na(MAF2))) {
            cat("Skipping row due to missing values.\n")
            next
          }

          ##=================Main-COLOC command================
          # Run coloc.abf
           my_res <- coloc.abf(
            dataset1 = list(snp = snp1, pvalues = p1, type = "quant", beta = Y1,  MAF = MAF1, N = N1),
            dataset2 = list(snp = snp2, pvalues = p2, type = "quant", beta = Y2, varbeta=V2, sdY=1)
          )
          ##===================================================

#          # Store results in the list
          result_row <- c(
            CHROM_y=CHROM_y, POS_y=POS_y, REF_y=REF_y, ALT_y=ALT_y, Index.Clumped=Index.Clumped, INDEX_ID=INDEX_ID, CHROM_x= CHROM_x, POS_x= POS_x, REF_x=REF_x, 
            ALT_x=ALT_x, rsid = rsid1, dataset2 = dataset2_file, coloc_result = my_res[[1]], geneSymbol = geneSymbol, tissue = tissue, se = se, nes = nes 
          )
          results_list[[length(results_list) + 1]] <- result_row
        }
      }
    }

    # Convert the results list to a dataframe
    results_df <- do.call(rbind, results_list)

    # Save the dataframe to CSV files with corresponding names in the output directory
    output_file <- file.path(output_directory, paste0(dataset_prefix, "_colocresults.csv"))
    write.table(results_df, file = output_file, row.names = FALSE,quote=FALSE)
  }
}
