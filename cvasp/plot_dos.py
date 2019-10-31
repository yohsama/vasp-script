#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import get_eigenvalue
import ase.io
def plot_dos(ISPIN):
    DOSCAR=open("DOSCAR")
    tmp=DOSCAR.readline()
    NIONS=int(tmp.split()[0])
    tmp=DOSCAR.readline()
    tmp=DOSCAR.readline()
    tmp=DOSCAR.readline()
    tmp=DOSCAR.readline()
    tmp=DOSCAR.readline()
    NEDOS=int(tmp.split()[2])
    TDOS=np.zeros((NEDOS,2*ISPIN+1))
    for i in range(NEDOS):
        tmp=DOSCAR.readline()
        TDOS[i]=np.array(tmp.split(),dtype=float)
    
    PDOS=[]
    for iNIONS in range(NIONS):
        tmp=DOSCAR.readline()    
        for i in range(NEDOS):
            tmp=DOSCAR.readline()
            PDOS.append(tmp.split())
    
    PDOS=np.array(PDOS,dtype=float)
    PDOS=PDOS.reshape((NIONS,NEDOS,PDOS.shape[1]))
    PDOS=PDOS[:,:,1:]
    
    EIG,OCC,FERMI=get_eigenvalue.get_eigenvalue()
    x=TDOS[:,0]-FERMI
    if ISPIN==2:
        plt.plot(x,np.array([TDOS[:,1],-TDOS[:,3]]).T,label="Total")
    else:
        plt.plot(x,TDOS[:,1],label="Total")
    atoms=ase.io.read("CONTCAR")
    symbol=atoms.get_chemical_symbols()
    for ele in set(symbol):
        select=[ele==i for i in symbol]
        if ISPIN==2:
            print(PDOS.shape)
            plt.plot(x,np.array([np.sum(PDOS[select,:,1::2],axis=(0,2)),-np.sum(PDOS[select,:,2::2],axis=(0,2))]).T,label=ele)
        else:
            print(PDOS.shape)
            plt.plot(x,np.sum(PDOS[select,:,1:],axis=(0,2)).T,label=ele)
    
    plt.xticks(fontsize=20)
    plt.xlim((-2,5))
    plt.legend()
    plt.show()

plot_dos(int(input('ISPIN\n')))
