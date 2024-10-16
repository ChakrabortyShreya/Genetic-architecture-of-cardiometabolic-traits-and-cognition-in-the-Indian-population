import sys
import os


path = sys.argv[1]
f = open(path,"r")
c =0
temp2=[]

for line in f:
	line = line.rstrip()
	if("---------------" not in line):
		temp = [i for i in line.split(" ") if i and "\n"!=line]
		if(temp!=[]):
			temp2.append(temp)
	
	else:
		index = 0
		if(len(temp2)<4):
			id = temp2[1][2].split(":")
			temp3 = "\t".join(["(INDEX)",temp2[1][2],"0","1.00",id[3],temp2[1][1],temp2[1][4]])
			#temp2.append(temp3)
			if("\n" in temp3):
				temp3 = temp3[:-1]
			print(temp3)

			temp2 = []
		else:
			for i in range(len(temp2)):
				if("(INDEX)" in temp2[i]):
					index = i
					temp2[i]="\t".join(temp2[i])
					if("\n" in temp2[i]):
						temp2[i] =temp2[i][:-1]
					print(temp2[i])

				if(index != 0 and i>index and "RANGE:" not in temp2[i]):
					temp2[i].insert(0," ")
					temp2[i]="\t".join(temp2[i])
					if("\n" in temp2[i]):
						temp2[i] =temp2[i][:-1]
					print(temp2[i])

				if("RANGE:" in temp2[i]):
					temp2=[]
					break
		print("-----------------------")

