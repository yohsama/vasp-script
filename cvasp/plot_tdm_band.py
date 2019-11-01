#!/usr/bin/python3
import ase.io
import numpy as np
import re
import matplotlib.pyplot as plt
from glob  import glob
from scipy.interpolate import griddata
import os
import sys
from multiprocessing import pool
from get_kpoints import *
from get_procar import get_procar
from get_wavecar import *
path='.'

def calc_tdm(coeff_1,coeff_2,eig_1,eig_2,igall_,b):
    #print('start')
    t1=time.time()
    dE=np.array(np.matrix(eig_1)-np.matrix(eig_2).T)/(2*13.605826)
    dE[dE==0]=np.inf
    tdm=np.zeros((dE.shape[0],dE.shape[1],4),dtype='complex')
    tdm_r=np.zeros((dE.shape[0],dE.shape[1],4),dtype='double')
    tdm_i=np.zeros((dE.shape[0],dE.shape[1],4),dtype='double')
    igall_=igall_.dot(b)
    for i in  range(3):
        #eigall_=np.zeros(igall_.shape)
        #eigall_[:,i]=igall_[:,i]
        #eigall_=eigall_.dot(b)
        #tdm_=(np.conj(coeff_1)*np.sum(eigall_,axis=1)).dot(coeff_2.T)
        tdm_=(np.conj(coeff_1)*(igall_[:,i])).dot(coeff_2.T)
        tdm_=1j/np.array(dE)*np.array(tdm_)*2.541746*0.529177249
#        tdm__r=np.real(tdm_)
#        tdm__i=np.imag(tdm_)
        if i==0:
            print("%f\t%f" % (np.real(tdm_),np.imag(tdm_)))
        tdm_= np.conj(tdm_)*tdm_
        tdm[:,:,i]=tdm_
#        tdm_r[:,:,i]=tdm__r
#        tdm_i[:,:,i]=tdm__i
    tmp=tdm[:,:,0]
    tdm[:,:,-1]=np.sum(np.array(tdm[:,:,0:3]),axis=2)
#    tdm_r[:,:,-1]=np.sum(np.array(tdm_r[:,:,0:3]),axis=2)
#    tdm_i[:,:,-1]=np.sum(np.array(tdm_i[:,:,0:3]),axis=2)
#    print(tdm_r)
#    print(tdm_i)
    #tdm=np.abs(tdm)
    #print('done')
    #print(time.time()-t1)
    return tdm





def cal_tdm_byband(ispin,iband,coeff,igall,eig,occ,b):
    tdm_k=np.array([])
    tdm_k_x=np.array([])
    tdm_k_y=np.array([])
    tdm_k_z=np.array([])
    for ikpt in range(coeff.shape[1]):
        coeff_=coeff[ispin,ikpt]
        eig_=eig[ispin,ikpt]
        igall_=igall[ispin,ikpt]
        tdm=calc_tdm(coeff_[iband[0]-1],coeff_[iband[1]-1],eig_[iband[0]-1],eig_[iband[1]-1],igall_,b)
        tdm_k=np.append(tdm_k,np.real(tdm[:,:,-1]))
        tdm_k_x=np.append(tdm_k_x,np.real(tdm[:,:,-4]))
        tdm_k_y=np.append(tdm_k_y,np.real(tdm[:,:,-3]))
        tdm_k_z=np.append(tdm_k_z,np.real(tdm[:,:,-2]))
    np.save('tdm_k.npy',tdm_k)
    #print(tdm_k_x)
    tdm=np.array([tdm_k_x,tdm_k_y,tdm_k_z,tdm_k]).reshape((4,coeff.shape[1])).T
    return tdm_k,tdm

def plot_tdm_band(kpoints,tdm_k,ispin,label=None,ax_tdm=None):
    #print(kpoints)
    ax_tdm.set_xlim([np.min(kpoints['kpath']),np.max(kpoints['kpath'])])
    ax_tdm.plot(kpoints['kpath'],tdm_k,label=r'$P_{('+str(label[0])+','+str(label[1])+')}$')
   # ax_tdm.plot(tdm_k,label=label)
    return ax_tdm


#from  plot_band import *
#
#kpoints=plot_band.from_kpoints('.',KPOINTS='IBZKPT',read=False,save=True)
#selectband=[[273,274]]
#ispin=0
#path='.'
#ax_tdm=plt.subplot()
#if True:
#    CONTCAR=path+'/CONTCAR'
#    WAVECARs=glob(path+'/WAVECAR_*')
#    if WAVECARs ==[] :
#        WAVECARs=[path+'/WAVECAR']
#    WAVECARs.sort()
#    print(WAVECARs)
#    tdm_k=np.array([])
#    for iselectband in selectband:
#        eig,occ=get_wavecar.get_eig(WAVECARs)
#        igall=get_wavecar.get_igall(WAVECARs,CONTCAR)
#        iband=np.arange(iselectband[0]-3,iselectband[1]+3)
#        coeff=get_wavecar.get_coeff(WAVECARs,iband)
#        print(coeff[:,:,iband])
#        print(igall.shape,coeff.shape)
#        tdm_k=cal_tdm_byband(ispin,iselectband,coeff,igall,eig,occ)
#        ax_tdm=plot_tdm_band(kpoints,tdm_k,0,label="_".join([str(i) for i in iselectband]),ax_tdm=ax_tdm)
#    print(tdm_k)
#    plt.legend()
#plt.tight_layout()
#plt.savefig('tband_%d.png' % ispin)


