import os
import pandas as pd
import sys

# argv[1]=Full path to *clumped_Formatted.txt
# argv[2]=Full path to .FamGrammarGamma.assoc_Processed
# argv[3]=Output path 

#Read data as a dataframe
data=pd.read_table(sys.argv[1], header=None,sep="\t")


#Read summary stats Processed file
sum_stat=pd.read_table(sys.argv[2],delim_whitespace=True)

#Store output file name
out=sys.argv[3]

#Remove lines starting with "------------"
data=data[data[0]!="-----------------------"]

#Will keep track of c=index and clumped SNPs
#SNP_ID_track = []

#To read and store chunks
#init_chunk=data.iloc[0].transpose() # Initialized dataframe which will store the current chunk
include_chunk=pd.DataFrame()
start_index=0
stop_string="(INDEX)"

#SNP_ID_track = pd.DataFrame()  

#Load chunks
def get_chunk(start_index):
	current_chunk=pd.DataFrame()
	counter=0
	while(start_index < len(data)):
		subset = data.iloc[start_index:start_index+1]
		if counter==0:
			current_chunk= pd.concat([current_chunk, subset])
			start_index+=1
			counter=1
			continue
		if subset[0].values != stop_string:
			current_chunk= pd.concat([current_chunk, subset])
			start_index+=1
		else:
			break
	
	return current_chunk,start_index

#To keep track on included SNPs
SNP_ID_track=[]

#Read chunks and include/exclude
while start_index<len(data):
	chunk=get_chunk(start_index) #Get chunk

	#if INDEX is not already clumped, then only include INDEX and their tagged variants
	if chunk[0][1].iloc[0]  not in SNP_ID_track:
		chunk[0]['Include']="Yes"
		SNP_ID_track.extend(chunk[0][1])
		
	else:
		chunk[0]['Include']="No"
	
	include_chunk=pd.concat([include_chunk,chunk[0]])
	start_index=chunk[1]

# Join summary stats on included set
Output=include_chunk[include_chunk['Include']=="Yes"]
Output.columns=["Index/Clumped","ID","KB","RSQ","ALLELES","F","P","Include"]

#Check for duplicates in included set
boolean = Output['ID'].duplicated().any()
if(boolean):
	print("Caution: Check for duplicates")
else:
	print("No worries :)")

new_Output=Output.merge(sum_stat[["ID","N_INFORMATIVE","AF","Beta","BetaVar","SE","Z"]], on="ID",how="left")
new_Output.to_csv(out,sep='\t', index=False)
