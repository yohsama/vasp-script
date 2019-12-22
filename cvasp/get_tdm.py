import numpy as np
import struct
import scipy
import scipy.integrate
import time
import ase.io
import matplotlib.pyplot as plt
from multiprocessing import Pool
from get_wavecar import get_wavecar 
path='./'
selectband=[[81,82]]

coeff,igall,eig,occ=get_wavecar('./WAVECAR',selectband)

def calc_tdm(coeff_,eig_,igall_):
    print('start')
    t1=time.time()
    #print(eig_)
    dE=np.array(np.matrix(eig_)-np.matrix(eig_).T)/(2*13.605826)
    dE[dE==0]=np.inf 
    tdm=np.zeros((4,dE.shape[0],dE.shape[1]),dtype='complex')
    for i in  range(3):
        tdm_=((np.conj(coeff_)*np.array(igall_)[:,i])).dot(coeff_.T)
        tdm_=1j/np.array(dE)*np.array(tdm_)*2.541746*0.529177249
        tdm_=tdm_*np.conj(tdm_)
        tdm[i]=tdm_
    tdm[-1]=np.sum(np.array(tdm),axis=0)
    tdm=np.abs(tdm)
    #print('done')
    #print(time.time()-t1)
    return tdm

spin=[0]
kpt=[0]
label=list(range(selectband[0][0],selectband[0][1]))
for ispin in spin:
    for ikpt in kpt:
        for iband in selectband:
            vb=int(np.sum(occ[ispin,ikpt])-iband[0])
            coeff_=coeff[ispin,ikpt]
            eig_=eig[ispin,ikpt]
            igall_=igall[ispin,ikpt]
            tdm=calc_tdm(coeff_[iband[0]-1:iband[1]],eig_[iband[0]-1:iband[1]],igall_)
      #      print(coeff_[iband[0]-1:iband[1]],eig_[iband[0]-1:iband[1]])
            np.save("tdm_"+str(ispin)+'_'+str(iband[0])+'_'+str(iband[1]),tdm)
       #     print(tdm)
            for id in range(4):
                plt.imshow(tdm[id])
                plt.yticks(np.arange(0,iband[1]-iband[0]),np.arange(iband[0],iband[1]))
                plt.xticks(np.arange(0,iband[1]-iband[0]),np.arange(iband[0],iband[1]))
                plt.show()
