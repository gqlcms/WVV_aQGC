import os
f1=open("filelistmc.txt")
dataset=[]
for i in f1:
	i=i.split(";")[0].replace(" ","").replace(".root","")
	dataset.append(i)
f2=open("dataset16.txt")
for j in f2:
	if j.replace("\n","") not in dataset:print j,"should be changed"
f2.close()
