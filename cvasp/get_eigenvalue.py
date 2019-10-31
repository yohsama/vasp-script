#!/usr/bin/env python
import numpy as np

def get_eigenvalue(file="EIGENVAL",OCC_E=0.001):
    EIGENVAL=open(file)
    tmp=EIGENVAL.readline()
    NIONS=int(tmp.split()[0])
    tmp=EIGENVAL.readline()
    tmp=EIGENVAL.readline()
    tmp=EIGENVAL.readline()
    tmp=EIGENVAL.readline()
    tmp=EIGENVAL.readline()
    NBAND=int(tmp.split()[2])
    NKPOINTS=int(tmp.split()[1])
    
    EIG=[]
    for ikpoints in range(NKPOINTS):
        tmp=EIGENVAL.readline()
        tmp=EIGENVAL.readline()
        for iband in range(NBAND):
            tmp=EIGENVAL.readline()
            EIG.append(tmp.split())
    
    EIG=np.array(EIG,dtype=float)
    EIG=EIG.reshape(NKPOINTS,NBAND,EIG.shape[1])
    if EIG.shape[2] == 3:
        OCC=EIG[:,:,2]
        EIG=EIG[:,:,1]
    else:
        OCC=EIG[:,:,3:]
        EIG=EIG[:,:,1:3]


    FERMI=np.max(EIG[OCC>OCC_E])
    return EIG,OCC,FERMI
