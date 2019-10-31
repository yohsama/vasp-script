#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import sys

LOC=open(sys.argv[1])
LOC_a=open(sys.argv[1]+"_alpha",'w')
LOC_b=open(sys.argv[1]+"_beta",'w')

tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
scale=float(tmp.split()[0])
cell=np.zeros((3,3))
for i in range(3):
    tmp=LOC.readline()
    LOC_a.write(tmp)
    LOC_b.write(tmp)
    cell[i,:]=tmp.split()[:3]

tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
NIONS=np.array(tmp.split(),dtype='int')
NIONS_TOT=np.sum(NIONS)
tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
for i in range(np.sum(NIONS)):
    tmp=LOC.readline()
    LOC_a.write(tmp)
    LOC_b.write(tmp)

tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
tmp=LOC.readline()
LOC_a.write(tmp)
LOC_b.write(tmp)
NGX,NGY,NGZ=np.array(tmp.split()[:3],dtype='int')
NG=(NGZ,NGY,NGX)
LOC_DAT=[]
for itmp in LOC:
    LOC_DAT.extend(itmp.split())

#LOC_DAT=np.array(LOC_DAT,dtype='double')
#LOC_DAT=[LOC_DAT.reshape(NG)]
LOC_DAT=[np.array(LOC_DAT[:np.prod(NG)],dtype='double'),np.array(LOC_DAT[-np.prod(NG):],dtype='double')]

LOC_DAT_a=(LOC_DAT[0]+LOC_DAT[1])/2
LOC_DAT_b=(LOC_DAT[0]-LOC_DAT[1])/2

count=0
for temp in LOC_DAT_a:
    count+=1
    print("%12.9f" % temp,end=" ",file=LOC_a)
    if(count%5==0):
        print("", file=LOC_a)

LOC_a.close()
for temp in LOC_DAT_b:
    count+=1
    print("%12.9f" % temp,end=" ",file=LOC_b)
    if(count%5==0):
        print("",file=LOC_b,flush=True)
LOC_b.close()


