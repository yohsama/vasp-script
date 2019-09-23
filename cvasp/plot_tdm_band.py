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
ISPIN=int(sys.argv[1])
setfermi=False
atoms=ase.io.read(path+'/CONTCAR')
pos=atoms.get_scaled_positions()
symbol=atoms.get_chemical_symbols()
rep=atoms.get_reciprocal_cell()

try:
    PRO=np.load(path+'/PRO.npy')
    EIG=np.load(path+'/EIG.npy')
    OCC=np.load(path+'/OCC.npy')
    KPOINTS=np.load(path+'/KPOINTS.npy')
    kpath=np.load(path+'/kpath.npy')
    kpt_label=np.load(path+'/kpt_label.npy')
except:
    i=0
    PROCARs=glob(path+'/PROCAR_*')
    if PROCARs == []:
        PROCARs=[path+"/PROCAR"]
    PROCARs.sort()
   # print(PROCARs)
    for PROCAR in PROCARs:
        EIG_,OCC_,PRO_,KPOINTS_=get_procar(PROCAR,ISPIN)
        if i==0:
            EIG=EIG_
            PRO=PRO_
            OCC=OCC_
            KPOINTS=KPOINTS_
        else:
            EIG=np.concatenate((EIG,EIG_),axis=1)
            PRO=np.concatenate((PRO,PRO_),axis=1)
            OCC=np.concatenate((OCC,OCC_),axis=1)
            KPOINTS=np.concatenate((KPOINTS,KPOINTS_),axis=0)
        i=i+1
    K_files=glob(path+'/KPOINTS_*')
    if K_files==[]:
        K_files=[path+"/KPOINTS"]
    K_files.sort()
    kpath,kpt_label=make_k_path(K_files,rep)
    np.save(path+'/PRO.npy',PRO)
    np.save(path+'/EIG.npy',EIG)
    np.save(path+'/OCC.npy',OCC)
    np.save(path+'/KPOINTS.npy',KPOINTS)
    np.save(path+'/kpt_label.npy',kpt_label)
    np.save(path+'/kpath.npy',kpath)


selectband=[[15,21]]

def calc_tdm(coeff_,eig_,igall_):
    #print('start')
    t1=time.time()
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

try:
    tdm_k=np.load(path+'/tdm_k.npy')
except:
    CONTCAR=path+'/CONTCAR'
    WAVECARs=glob(path+'/WAVECAR_*')
    if WAVECARs ==[] :
        WAVECARs=[path+'/WAVECAR']
    WAVECARs.sort()
    tdm_k=[[],[]]
    for WAVECAR in WAVECARs:
        coeff,igall,eig,occ=get_wavecar(WAVECAR,CONTCAR,selectband)
        for ispin in range(ISPIN):
            for ikpt in range(coeff.shape[1]):
                #print(coeff.shape)
                for iband in selectband:
                    coeff_=coeff[ispin,ikpt]
                    eig_=eig[ispin,ikpt]
                    igall_=igall[ispin,ikpt]
                   # print(iband)
                    tdm=calc_tdm(coeff_[iband[0]-1:iband[1]],eig_[iband[0]-1:iband[1]],igall_)
                tdm_k[ispin].append(tdm[-1,0,-1])
    np.save('tdm_k.npy',tdm_k)

if setfermi==True:
    EIG_=EIG
    EIG_[OCC<0.01]=0
    Nb=np.prod(EIG_.shape[1:])
    fermi=np.max(EIG_)
    fermi_ax=np.argmax(EIG_.reshape((ISPIN,Nb)))
    print("Fermi in spin %s is %f" % (("alpha","beta")[int(fermi_ax//Nb)],fermi))
    EIG=EIG-fermi

EIG=EIG-7.9

L_INX=[(i=="Si") for i in symbol]
if kpath.shape[0]==1:
    PRO=np.concatenate((PRO,PRO),axis=1)
    EIG=np.concatenate((EIG,EIG),axis=1)
    kpath=np.array([0,1])
    intx=np.linspace(0,np.max(kpath),101)
    fig=plt.figure(figsize=(3*ISPIN,10))
else:
    intx=np.linspace(0,np.max(kpath),1001)
    fig=plt.figure(figsize=(8*ISPIN,10))

Elim=(-5,10)
EIG_=np.max(np.max(EIG,axis=0),axis=0)
RL=np.argmax(EIG_>Elim[0])
EIG_=np.min(np.min(EIG,axis=0),axis=0)
RH=np.argmax(EIG_>Elim[1])
if RH==0:
    RH=EIG.shape[2]
print(RL,RH)
PRO=PRO[:,:,RL:RH]
EIG=EIG[:,:,RL:RH]




from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec


for ispin in range(ISPIN):
    print(np.sum(PRO[ispin],axis=2).shape)
    L1=np.sum(PRO[ispin,:,:,L_INX],axis=0)/np.sum(PRO[ispin],axis=2)
    L=L1#/(L1+L2)
    c=griddata(kpath,np.array(L),intx)
    E=griddata(kpath,EIG[ispin,:],intx)
    print((np.array([intx.T]*(EIG.shape[1])).T.shape,c.shape,E.shape))
    if ISPIN==1:
        ax=plt.subplot(111)
        divider=make_axes_locatable(ax)
        ax_cbar = divider.append_axes('right', size="5%", pad=0.02)
        sc=ax.scatter(np.array([intx.T]*(EIG.shape[2])).T,E,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='jet')
        if any(kpt_label!=""):
            ax.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
        #    ax.set_xticks(kpath[kpt_label!=""])
        #    ax.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
        else:
            ax.set_xticks([])
        ax.set_xlim([np.min(kpath),np.max(kpath)])
        ax.set_ylim([Elim[0],Elim[1]])
        ax.set_yticklabels(ax.get_yticks(),fontsize=20)
        cb=plt.colorbar(sc,cax=ax_cbar,orientation='vertical')
        cb.set_ticks([0,1])
        cb.ax.tick_params(labelsize=20)
        cb.set_ticklabels(['C','Si'])
        ax_tdm= divider.append_axes('bottom', size="20%", pad=0.1)
        if any(kpt_label!=""):
            ax_tdm.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
            ax_tdm.set_xticks(kpath[kpt_label!=""])
            ax_tdm.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
        ax_tdm.set_xlim([np.min(kpath),np.max(kpath)])
        ax_tdm.set_ylim([-2,np.max(tdm_k[ispin])])
        ax_tdm.plot(kpath,tdm_k[ispin])
    else:
        gs=GridSpec(1,9)
        if ispin==1:
            ax1=plt.subplot(gs[0,0:4])
            sc1=ax1.scatter(np.array([intx.T]*(EIG.shape[2])).T,E,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='jet')
            if any(kpt_label!=""):
                ax1.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
                ax1.set_xticks(kpath[kpt_label!=""])
                ax1.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
            else:
                ax1.set_xticks([])
            ax1.set_xlim([np.min(kpath),np.max(kpath)])
            ax1.set_ylim([-2,5])
            ax1.set_yticklabels(ax1.get_yticks(),fontsize=20)
        else:
            ax2 = plt.subplot(gs[0,4:8])
            #ax2.set_yticks([])
            divider2=make_axes_locatable(ax2)
            #ax_cbar=divider2.append_axes('right', size="5%", pad=0.2)
            sc2=ax2.scatter(np.array([intx.T]*(EIG.shape[2])).T,E,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='jet')
            ax2.set_yticks([])
            if any(kpt_label!=""):
                ax2.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
                ax2.set_xticks(kpath[kpt_label!=""])
                ax2.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
            else:
                ax2.set_xticks([])
            ax2.set_xlim([np.min(kpath),np.max(kpath)])
            ax2.set_ylim([-2,5])
            ax_cbar=plt.subplot(gs[:,8])
            cb=plt.colorbar(sc2,cax=ax_cbar,orientation='vertical')
            cb.set_ticks([0,1])
            cb.ax.tick_params(labelsize=20)
            cb.set_ticklabels(['C','Si'])
        plt.subplots_adjust(wspace=0.02)
plt.tight_layout()
plt.savefig('band.png',dpi=360)
plt.show(cb)

