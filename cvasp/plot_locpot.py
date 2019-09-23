#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import sys

LOC=open('LOCPOT')

tmp=LOC.readline()
tmp=LOC.readline()
scale=float(tmp.split()[0])
cell=np.zeros((3,3))
for i in range(3):
    tmp=LOC.readline()
    cell[i,:]=tmp.split()[:3]

tmp=LOC.readline()
tmp=LOC.readline()
NIONS=np.array(tmp.split(),dtype='int')
NIONS_TOT=np.sum(NIONS)
tmp=LOC.readline()
for i in range(np.sum(NIONS)):
    tmp=LOC.readline()

tmp=LOC.readline()
tmp=LOC.readline()
NGX,NGY,NGZ=np.array(tmp.split()[:3],dtype='int')
NG=(NGZ,NGY,NGX)
LOC_DAT=[]
for itmp in LOC:
    LOC_DAT.extend(itmp.split())

try:
    LOC_DAT=np.array(LOC_DAT,dtype='double')
    LOC_DAT=[LOC_DAT.reshape(NG)]
except:
    LOC_DAT=[np.array(LOC_DAT[:np.prod(NG)],dtype='double').reshape(NG),np.array(LOC_DAT[-np.prod(NG):],dtype='double').reshape(NG)]

_axe=['c','b','a'].index(sys.argv[1])
#_axe=int(sys.argv[1])-1
axe=[0,1,2]
axe.remove(_axe)
axe=tuple(axe)
LOC_DAT=np.array(LOC_DAT)
print(LOC_DAT.shape)
if LOC_DAT.shape[0]>1:
    LOC_DAT_A=(LOC_DAT[0]+LOC_DAT[1])/2
    LOC_DAT_B=(LOC_DAT[0]-LOC_DAT[1])/2
else:
    LOC_DAT_A=(LOC_DAT[0])
plt.plot(range(NG[_axe])/NG[_axe]*np.linalg.norm(cell[(2,1,0)[_axe],:]),np.average(LOC_DAT_A,axis=axe))
#try:
#    plt.plot(range(NG[_axe])/NG[_axe]*np.linalg.norm(cell[(2,1,0)[_axe],:]),np.average(LOC_DAT_B,axis=axe))
#except:
#    pass

plt.savefig("LOC_%s.jpg" % sys.argv[1])
plt.show()

