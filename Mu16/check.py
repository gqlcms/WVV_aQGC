import os
file=os.listdir("/eos/cms/store/user/xulyu/data2016_9_8_corr_ptwwong")

print file,len(file)
iin=[]

f=open("filelist.txt")
for l in f:
	l=l.split(";")
	l=l[0].replace(" ","")
	iin.append(l);
for i in file:
	if i not in iin:print i
