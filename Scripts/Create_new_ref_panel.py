import os
from mpi4py import MPI

#chromosomes = [f'{i}' for i in range(1, 23)]
#chrom_num=[f'{i}' for i in range(1, 23)]

thread_count = str(os.environ.get('SLURM_CPUS_PER_TASK', 1))
ref_path="/Ref-Panel-TLSA-SAS"

def sas_tlsa_panel(num):
	chrom_dir = os.path.join(os.getcwd(), ref_path, chromosomes[num])
  
  # Chunk chromomomes to merge TLSA and SAS WGS in chunks
	with open(f"Chromosome_Lengths.txt","r")as lengths:
		for line in lengths:
			cols=line.strip().split('\t')
			if cols[0]==chromosomes[num]:
				chr_len=int(cols[1])
	chunk_size= 5000000 #divide chromosome into chunk size of 5KB
	num_chunks_list = []
	start_list = []
	end_list = []
	# Calculate the number of chunks
	num_chunks = chr_len // chunk_size
	if chr_len % chunk_size != 0:
		num_chunks += 1
	num_chunks_list = list(range(1, num_chunks + 1))
	# Calculate the starting and ending positions of each chunk
	for i in range(0, chr_len, chunk_size):
		start = i+1
		end = min(i + chunk_size, chr_len)
		start_list.append(start)
		end_list.append(end)
		
	os.system(f"time mkdir -p {chrom_dir}/out")
	os.system(f"time mkdir {chrom_dir}/log")
	with open(f"{chrom_dir}/{chromosomes[num]}-Chunks.txt","w") as file:
		file.write(f"Chromosome Number: {chromosomes[num]}\n")
		file.write(f"Chromosome Length: {chr_len}\n\n")
		file.write("Chunk ID\tStart BP\tEnd BP\n")
		for chunk_id, start_bp, end_bp in zip(num_chunks_list, start_list, end_list):
			file.write(f"{chunk_id}\t\t{start_bp}\t\t{end_bp}\n")

	#Converting the QC-d TLSA-vcf.gz datasets into .hap/.legend formats
	os.system(f"time bcftools convert --threads 72 -h {chromosomes[num]}/{chromosomes[num]} {chrom_dir}/{chromosomes[num]}_QCd.phased_SNVs_INDELs_sorted.vcf.gz")
	#Unzipping the compressed hap legend files
	os.system(f"time bgzip -d {chrom_dir}/{chromosomes[num]}.hap.gz")
	os.system(f"time bgzip -d {chrom_dir}/{chromosomes[num]}.legend.gz")

	#TO PRODUCE SAS-hap and legend files from new QC'd SAS WGS files:
  os.system(f"time bcftools convert --threads 72 -h {sas_dir}/{chromosomes[num]}/tmp-SAS-{chromosomes[num]} {sas_dir}{chrom_dir}/SAS_{chromosomes[num]}_QCd.phased_SNVs_INDELs_sorted.vcf.gz")
	#Unzipping the compressed hap legend files
	os.system(f"time bgzip -d {sas_dir}/{chrom_dir}/tmp-SAS-{chromosomes[num]}.hap.gz")
	os.system(f"time bgzip -d {sas_dir}/{chrom_dir}/tmp-SAS-{chromosomes[num]}.legend.gz")
	'''	

	#Making output and log files for chunk generation and logging respectively
	#Creating hap/legend chunks upon merging SAS and TLSA files
	for i in range(num_chunks):
		#For each chromosome demarcated by chromosomes[num] , the command would create as many chunks as i, with start_list[i] and end_list[i] depicting the boundaries of the region to be merged.
		print("Chunk ",num_chunks_list[i]," for ",chromosomes[num],"starting at ",start_list[i],"ending at ",end_list[i])
		#Creating hap/legend chunks upon merging SAS and TLSA files
		os.system(f"time impute2 -m {map_dir}/{chromosomes[num]}.b38.gmap.gz -h {sas_ref_dir}/{chromosomes[num]}/tmp-sas-{chromosomes[num]}.hap  {chrom_dir}/{chromosomes[num]}.hap  -l {sas_ref_dir}/{chromosomes[num]}/tmp-sas-{chromosomes[num]}.legend {chrom_dir}/{chromosomes[num]}.legend -merge_ref_panels_output_ref {chrom_dir}/out/{chromosomes[num]}-chunk-{str(num_chunks_list[i])}  -r {chrom_dir}/log/chunk_{str(num_chunks_list[i])}.log -int {str(start_list[i])} {str(end_list[i])}")
	
#Merging the reference panels chunk-wise
	os.system(f"for i in $(ls -1v {chrom_dir}/out/{chromosomes[num]}-chunk-*.legend);do cat $i >> {chrom_dir}/out/tmp-{chromosomes[num]}-merge-chunks.legend ;done")
	os.system(f"for i in $(ls -1v {chrom_dir}/out/{chromosomes[num]}-chunk-*.hap);do cat $i >>{chrom_dir}/out/{chromosomes[num]}-merge-chunks.hap;done")
	with open(f"{chrom_dir}/out/tmp-{chromosomes[num]}-merge-chunks.legend","r") as f_in ,open(f"{chrom_dir}/out/{chromosomes[num]}-merge-chunks.legend","w") as f_out :
		header=f_in.readline()
		f_out.write(header)
		# remove all other occurrences of the header and write the remaining content to the output .legend file
		for line in f_in:
			if line != header:
				f_out.write(line)
	os.system(f"time bgzip -@ 72 {chromosomes[num]}/out/{chromosomes[num]}-merge-chunks.hap ")
	os.system(f"time bgzip -@ 72 {chromosomes[num]}/out/{chromosomes[num]}-merge-chunks.legend")
	os.system(f"rm {chromosomes[num]}/out/tmp-*")
	os.system(f"mv {chromosomes[num]}/out/{chromosomes[num]}-merge-chunks* {chromosomes[num]}/")
	os.system(f"bcftools query -l {chrom_dir}/{chromosomes[num]}_QCd.phased_SNVs_INDELs_sorted.vcf.gz > {chromosomes[num]}/tmp-{chromosomes[num]}-TLSA.samples")
	with open(f"{sas_ref_dir}/{chromosomes[num]}/{chromosomes[num]}-SAS.samples",'r') as f1, open(f"{chromosomes[num]}/tmp-{chromosomes[num]}-TLSA.samples",'r') as f2, open(f"{chromosomes[num]}/{chromosomes[num]}-merge-chunks.samples",'w') as f_out:
		# read sample IDs from file1
		ids1 = f1.read().splitlines()
		# read sample IDs from file2
		ids2 = f2.read().splitlines()
		# paste the IDs twice side by side as columns separated by a space
		id_pairs1 = [f"{id1} {id1}" for id1 in ids1]
		id_pairs2 = [f"{id2} {id2}" for id2 in ids2]
		# concatenate the two new files and write the ID header line and the concatenated ID pairs to the output file
		f_out.write('ID1 ID2\n')
		f_out.write('\n'.join(id_pairs1 + id_pairs2))
  
#Converting hap/legend files into vcf.gz
	os.system(f"bcftools convert --threads 72 -H {chromosomes[num]}/{chromosomes[num]}-merge-chunks -Oz -o {chromosomes[num]}/{chromosomes[num]}-reference_panel.vcf.gz")
	#os.system(f"rm {chromosomes[num]}/tmp-*")
	
  os.system(f"time bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' {chromosomes[num]}/{chromosomes[num]}-reference_panel.vcf.gz  --threads 72 -Oz -o  {chromosomes[num]}/{chromosomes[num]}-reference_panel-with-IDs.vcf.gz")
	os.system(f"bcftools index -t --threads 72 {chromosomes[num]}/{chromosomes[num]}-reference_panel-with-IDs.vcf.gz")
	os.system(f"bcftools sort {chromosomes[num]}/{chromosomes[num]}-reference_panel-with-IDs.vcf.gz -Oz -o {chromosomes[num]}/{chromosomes[num]}-mac2-merged-sorted-reference_panel-with-IDs.vcf.gz") 
	os.system(f"bcftools index -t --threads 72 {chromosomes[num]}/{chromosomes[num]}-mac2-merged-sorted-reference_panel-with-IDs.vcf.gz")

  

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


if rank < (len(chromosomes)):
	#sas_tlsa_panel(rank)
