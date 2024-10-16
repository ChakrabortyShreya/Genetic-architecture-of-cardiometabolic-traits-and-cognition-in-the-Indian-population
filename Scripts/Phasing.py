# Define the list of chromosomes using a more readable format
chromosomes = [f'{i}' for i in range(1, 23)]
chrom_num=[f'{i}' for i in range(1, 23)]

base_name="phase_input_without extension"
mapdir='chr-maps'
thread_count = str(os.environ.get('SLURM_CPUS_PER_TASK', 1))

def phasing(num):
	chrom_dir = os.path.join(os.getcwd(), chromosomes[num])
	if not os.path.exists(chrom_dir):
		os.makedirs(chrom_dir)
	begin_time = start_time = time.time()
	log_path = os.path.join(chrom_dir, f"{chromosomes[num]}_ShapeIt5_RunTime.log")
	log_file=open(log_path, 'w')
	log_file.write(f"RUNTIME REPORT - {chromosomes[num]}\n\nAnnotation\n")
	#Common variant phasing (MAF > 0.001 [0.1%])
	main_file=f'{chrom_dir}/{chromosomes[num]}_{base_name}.vcf.gz'
	os.system(f"phase_common --input {main_file}  --filter-maf 0.001 --region {chrom_num[num]}  --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --output {chrom_dir}/{chromosomes[num]}_{base_name}_phased.vcf  --thread {thread_count}")
	os.system(f"bgzip -@ {thread_count} {chrom_dir}/{chromosomes[num]}_{base_name}_phased.vcf")
	os.system(f"bcftools index -t --threads {thread_count} {chrom_dir}/{chromosomes[num]}_{base_name}_phased.vcf.gz")
	end_time = time.time()
	elapsed_time = round((end_time - start_time) / 60, 3)
	log_file.write(f"\tTook {elapsed_time} minutes\n")

def rare_phasing(num):
	chrom_dir = os.path.join(os.getcwd(), chromosomes[num])
	if not os.path.exists(chrom_dir):
		os.makedirs(chrom_dir)
	main_file=f'{chrom_dir}/{chromosomes[num]}_{base_name}.vcf.gz'
	os.system(f"GLIMPSE2_chunk --input {main_file} --region {chrom_num[num]}  --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --sequential --output {chrom_dir}/chunks_{chromosomes[num]}.txt")
	os.system(f"mkdir {chrom_dir}/chunks")
	with open(f'{chrom_dir}/chunks_{chromosomes[num]}.txt', 'r') as file:
		for line in file: 
			os.system(f'phase_rare --input {main_file} --scaffold {chrom_dir}/{chromosomes[num]}_{base_name}_phased.vcf.gz --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --input-region {line.split()[3]} --scaffold-region {line.split()[2]} --output {chrom_dir}/chunks/{chromosomes[num]}_{base_name}.target.phased.chunk{line.split()[0]}.vcf --thread {thread_count}')
	os.system(f"ls -1v  {chrom_dir}/chunks/{chromosomes[num]}_{base_name}.target.phased.chunk*.vcf >  {chrom_dir}/chunks/{chromosomes[num]}_chunks_list.txt")
	os.system(f"bcftools concat -Oz --threads {thread_count} -o {chrom_dir}/{chromosomes[num]}_{base_name}.target.phased.vcf.gz -f {chrom_dir}/chunks/{chromosomes[num]}_chunks_list.txt")
	os.system(f"bcftools index -t --threads {thread_count}  {chrom_dir}/{chromosomes[num]}_{base_name}.target.phased.vcf.gz")



comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank < len(chromosomes):
	phasing(rank)
	rare_phasing(rank)
