#coding:utf-8

import sys
import os

redseq=open(sys.argv[1],"r")
rmk=open(sys.argv[1]+".tmp","w")
num=int(sys.argv[2])
n=0

for i in redseq:
	i = i.strip()
	if n <= 2:
		n+=1
		rmk.write(i+"@@KIREME@@")
	elif n == 3:
		n = 0
		rmk.write(i+"\n")
redseq.close()
rmk.close()

tmp=open(sys.argv[1]+".tmp","r")

for t in tmp:
	t = t.strip()
	tt = t.split("@@KIREME@@")
	seq=tt[1]
	if len(seq) >= num:
		print(tt[0]+"\n"+tt[1]+"\n"+tt[2]+"\n"+tt[3])
tmp.close()

os.system("rm "+sys.argv[1]+".tmp")
