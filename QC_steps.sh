#!/bin/bash

read -p 'Enter vcf filename with path ("../1.0.Merged/SANSCOG_2_108_v1.vcf.gz" kept in the merge directory):' file
read -p ' Sample list to be excluded (Exclude_Samples.txt):' exclude
read -p 'Enter output file name ("SANSCOG_2_108_QC_d"):' output 
read -p 'Enter output name for GRCh37 format ("SANSCOG_2_108_QC_d_chr_renamed"):' renamed

mkdir Temp

module load plink/2
module load R/4.0.5

##Extract Autosomes and biallelic variants
plink2 --vcf $file --remove Required/$exclude --chr 1-22,X,Y,XY,MT,PAR1,PAR2 --max-alleles 2 --min-alleles 2  --make-pgen --out Temp/Intermediate_file1 --allow-extra-chr 

## Sample missingness - > 10% --remove
plink2 --pfile Temp/Intermediate_file1 --make-pgen --out Temp/Intermediate_file1.1 --mind 0.1


##Impute and remove ambiguous ID
module load plink/1.9.0
plink2 --pfile Temp/Intermediate_file1.1 --make-bed --out Temp/Intermediate_file1.1.1

plink --bfile Temp/Intermediate_file1.1.1 --impute-sex --out Temp/Intermediate_file2 --make-bed
plink --bfile Temp/Intermediate_file2 --check-sex  --out Temp/Intermediate_file2.0


#Keep ambiguous
plink2 --pgen Temp/Intermediate_file1.1.pgen --pvar Temp/Intermediate_file1.1.pvar  --psam Temp/Intermediate_file2.fam --out Temp/Intermediate_file1.2 --make-pgen --chr 1-22


##Remove duplicate variants
#plink2 --pfile Temp/Intermediate_file1.2 --set-all-var-ids @:#[b38]\$r,\$a --rm-dup force-first --make-pgen --out Temp/Intermediate_file4 --new-id-max-allele-len 2516 
plink2 --pfile Temp/Intermediate_file1.2 --rm-dup force-first --make-pgen --out Temp/Intermediate_file4 --new-id-max-allele-len 2516



plink2 --pfile Temp/Intermediate_file4 --geno 0.02 --out Temp/Intermediate_file5 --make-pgen

#plink2 --pfile Temp/Intermediate_file5 --mac 3 --out Temp/Intermediate_file6 --make-pgen

#plink2 --pfile Temp/Intermediate_file6 --hwe 1e-7 --out QC_d_file --make-pgen

grep -v "#" Temp/Intermediate_file5.pvar |cut  -f 3 > Temp/VarList_Int_QC_d_in.txt

grep -v "#" Temp/Intermediate_file5.psam |cut  -f 1 > Samples_retained_after_QC.txt

##Extract samples and variants retained after QC
module load bcftools-1.11.3

#time bcftools annotate --set-id '%CHROM\:%POS\[b38\]%REF\,%ALT' -Oz -o Temp/inVCF.id.renamed.vcf.gz --threads 72 $file 
#time bcftools view --include ID==@VarList_QC_d_in.txt -S Samples_retained_after_QC.txt Temp/inVCF.id.renamed.vcf.gz -Oz -o $output.vcf.gz --threads 72

time bcftools view --include ID==@Temp/VarList_Int_QC_d_in.txt -S Samples_retained_after_QC.txt $file -Oz -o Temp/Intermediate_file5.vcf.gz --threads 72
time bcftools index -t Temp/Intermediate_file5.vcf.gz --threads 72


#MAF and  Hardy fill tags >=3
time bcftools +fill-tags Temp/Intermediate_file5.vcf.gz -Oz -o Temp/Intermediate_file6.vcf.gz --threads 72 -- -t HWE,F_MISSING,MAF
time bcftools index -t Temp/Intermediate_file6.vcf.gz --threads 72

time bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MAF\t%INFO/AN\t%INFO/HWE\n' Temp/Intermediate_file6.vcf.gz > Temp/MAF_AN_Hardy_Extract.txt
time python Required/MAC_Hardy_Extract.py


time bcftools view --include ID==@Temp/MAC2_Hardy_in.txt Temp/Intermediate_file6.vcf.gz -Oz -o $output.MAC2.vcf.gz --threads 72

time bcftools index -f -t $output.MAC2.vcf.gz --threads 72



##Rename chromosomes
time bcftools annotate $output.MAC2.vcf.gz --rename-chrs Required/rename_chr.txt -Oz -o $renamed.MAC2.vcf.gz --threads 72
time bcftools index -f -t $renamed.MAC2.vcf.gz 
