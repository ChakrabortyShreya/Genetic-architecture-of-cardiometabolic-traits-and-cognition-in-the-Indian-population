import subprocess
import gzip
import os
from mpi4py import MPI

#performing for Chr1-22
#x=[f'{i}' for i in range(1,23)]
#y=[str(i) for i in range(1,23)]


base_dir=''
ref_dir=f"{base_dir}/Phasing/"   ## Where reference panel has been created
map_dir = f"{base_dir}/chr-maps"            ## Genetic maps 
sanscog_dir="SANSCOG/mac2"  #path to genotyped data
imp_dir=f"{base_dir}/Imputation_Result"

A=str("A")
T=str("T")
C=str("C")
G=str("G")

thread_count = str(os.environ.get('SLURM_CPUS_PER_TASK', 1))


def imputation(num):
	os.system(f"mkdir -p {x[num]}/out")
	os.system(f"mkdir {x[num]}/log")
	#Rename chromosomes
  os.system(f"bcftools annotate {sanscog_dir}/chr{x[num]}/SANSCOG-chr{x[num]}-QC-phased-wIDs-mac2-annot.vcf.gz --rename-chrs chr_map.txt -Oz -o SANSCOG/SANSCOG-chr{x[num]}-QC-phased-wIDs-mac2-annot_chrrenamed.vcf.gz --threads {thread_count}")
	os.system(f"bcftools index -t SANSCOG/SANSCOG-chr{x[num]}-QC-phased-wIDs-mac2-annot_chrrenamed.vcf.gz --threads {thread_count}")
  
  #Chunk chromosomes for imputation
 	os.system(f"imp5Chunker_1.1.5_static --h {ref_dir}/{x[num]}/chr{x[num]}_reference_panel.target.phased.vcf.gz --g SANSCOG/SANSCOG-chr{x[num]}-QC-phased-wIDs-mac2-annot_chrrenamed.vcf.gz  --window-size 5000000  --r chr{y[num]}  --o {x[num]}/{x[num]}-coordinates.txt")
	buffer_region=[]
	imp_region=[]
	with open(f"{x[num]}/{x[num]}-coordinates.txt",'r') as chunks:
		for line in chunks:
			cols=line.strip().split('\t')
			#Saving 3rd column Buffer Region
			buffer_region.append(cols[2])
			#Saving 4th column Imputation Region
			imp_region.append(cols[3])
	print("Buffer Region\n",buffer_region,"\n Imputation Regions\n ",imp_region)   
	
  ## Imputation process split across multiple chunks according to the Chromosomal Map
	c=0
	while (c<len(buffer_region)):
		z=str(c+1)
		os.system(f"time impute5_1.1.5_static --m {map_dir}/chr{x[num]}.b38.gmap.gz --h {ref_dir}/{x[num]}/chr{x[num]}_reference_panel.target.phased.vcf.gz --g SANSCOG/SANSCOG-chr{x[num]}-QC-phased-wIDs-mac2-annot_chrrenamed.vcf.gz  --r {imp_region[c]}  --o {x[num]}/out/{x[num]}-chunk-{z}-imputed.vcf.gz --l {x[num]}/log/{x[num]}-chunk-{z}.log --threads {thread_count} --out-gp-field --buffer-region {buffer_region[c]}")
		c+=1
	## Merging all the chunks in order
	os.system(f"ls -1v {x[num]}/out/{x[num]}-chunk-*-imputed.vcf.gz > {x[num]}/{x[num]}-chunks_list.txt")
	## Concatenation
	os.system(f"bcftools concat {x[num]}/out/{x[num]}-chunk-*-imputed.vcf.gz -Oz -o {x[num]}/SANSCOG-{x[num]}-imputed-mac2.vcf.gz --threads {thread_count}")
	## Sorting
	os.system(f"bcftools sort {x[num]}/SANSCOG-{x[num]}-imputed-mac2.vcf.gz -Oz -o {x[num]}/SANSCOG-{x[num]}-imputed-sorted-mac2.vcf.gz")
	os.system(f"bcftools index -t {x[num]}/SANSCOG-{x[num]}-imputed-sorted-mac2.vcf.gz --threads {thread_count}")
	## Adding MAF in INFO field
	os.system(f"bcftools +fill-tags --threads {thread_count} {x[num]}/SANSCOG-{x[num]}-imputed-sorted-mac2.vcf.gz  -Oz -o {x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz -- -t MAF")
	os.system(f"bcftools index -t  {x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz --threads {thread_count}")
	## Deleting the intermediate chunk files
	#os.system(f"rm -rf {x[num]}/out/")
	## Extracting the Variant ID, Allele Frequency and R.Sq. Values for Metric File creation
	os.system(f"bcftools query -f '%ID\t%INFO/AF\t%INFO/INFO\n' {x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz > {x[num]}/SANSCOG-{x[num]}-unfiltered-Metrics.txt")
	os.system(f"bgzip -@ {thread_count} {x[num]}/SANSCOG-2_108-{x[num]}-unfiltered-Metrics.txt")

	### The following commented line is to produce a filtered file with R.Sq > 0.6 and ATGC Flip filter
	os.system(f"bcftools view -e '(INFO/INFO <0.6 & INFO/MAF>=0.4 & INFO/MAF<=0.6 & ((REF=\"A\" & ALT=\"T\") | (REF=\"T\" & ALT=\"A\")| (REF=\"G\" & ALT=\"C\")| (REF=\"C\" & ALT=\"G\")))'  --threads {thread_count} {x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz -Oz -o {x[num]}/SANSCOG-{x[num]}-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem.vcf.gz")
	os.system(f"bcftools index -t --threads {thread_count} {x[num]}/SANSCOG-{x[num]}-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem.vcf.gz")

##Postimputation processing for plots
def postimp(num):
	main_file=f'{x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz'
	##Getting sample count for MAC Calculation
	command=f"bcftools query -l {main_file} |wc -l"
	sample_num= subprocess.run([command],text=True,stdout=subprocess.PIPE,shell=True,stderr=subprocess.PIPE)
	sample_count=float(sample_num.stdout)
#	print(sample_count)
	## Calculating MAC and MAF using python built in commands
	with gzip.open(f"{x[num]}/SANSCOG-{x[num]}-unfiltered-Metrics.txt.gz",'rt') as ip, gzip.open(f"{x[num]}/SANSCOG-{x[num]}-unfiltered-MAF-MAC-Metrics.txt.gz",'wt') as op_mac:
		for line in ip:
			cols=line.strip().split('\t')
			if float(cols[1]) > 0.5:
				maf_new=float(1-float(cols[1]))
				cols[1]=str(maf_new)
			mac_new = round((float(cols[1]) * (2*sample_count)))
			output = f"{cols[0]}\t{cols[1]}\t{mac_new:.0f}\t{cols[2]}\n"
			op_mac.write(output)
	metrics_file= gzip.open(f'{x[num]}/SANSCOG-{x[num]}-unfiltered-MAF-MAC-Metrics.txt.gz','rt')
	op6maf2ids= gzip.open(f'{x[num]}/SANSCOG-{x[num]}-info0.6-flip-0.0002_IDs.txt.gz','wt')
	next(metrics_file)
	for line in metrics_file:
		cols=line.strip().split('\t')
		id=cols[0].strip().split(':')
		ref=id[2]
		alt=id[3]
		maf=float(cols[1])
		info=float(cols[3])
		##The following is the criteria for Allele Flip(ATGC)
		flip_check = (0.4 <= maf <= 0.6) and ((ref in ['A','T'] and alt in ['A','T']) or (ref in ['G','C'] and alt in ['G','C']))
		if (info >= 0.6 and maf>=0.0002 and not flip_check):
			op6maf2ids.write(str(id[0])+'\t'+str(id[1])+'\n')

##Extract variants with imputation accuracy >0.6, no A/T G/C snps with AF between 40-60% and MAF>0.02% 
def extraction(num):
	print(f"{x[num]}/SANSCOG-{x[num]}_info0.6-flip-0.0002_IDs.txt.gz")
	os.system(f"bcftools view -R  {x[num]}/SANSCOG-{x[num]}--info0.6-flip-0.0002_IDs.txt.gz {x[num]}/SANSCOG-{x[num]}-imputed-unfiltered-sorted-mac2-annot.vcf.gz --threads {thread_count} -Oz -o {x[num]}/SANSCOG-{x[num]}-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-extract.vcf.gz")
	os.system(f"bcftools index -t --threads {thread_count} {x[num]}/SANSCOG{x[num]}-imputed-filtered-0.6-sorted-mac2-annot-atgc-rem-extract.vcf.gz")


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank < (len(x)):
	imputation(rank)
	postimp(rank)
	extraction(rank)

