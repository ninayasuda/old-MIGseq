#coding:utf-8
#作成日：2018/OCT/02
#作成者：湯淺英知(東工大・生命理工)
#使い方：python make_list_1st.py <mappingのidのlist> > out.list

import sys

fil=open(sys.argv[1],"r")
box=[]

for i in fil:
	i = i.strip()
	if len(i) > 0:
		ii = i.split("_Sample_")
		id = ii[0]
		if id not in box:
			print(id)
			box.append(id)
fil.close()

