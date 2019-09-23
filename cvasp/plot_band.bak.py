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

path='.'
ISPIN=int(sys.argv[1])
setfermi=True
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
    print(PROCARs)
    for PROCAR in PROCARs:
        EIG_,OCC_,PRO_,KPOINTS_=get_procar(PROCAR,ISPIN)
        if i==0:
            EIG=EIG_
            PRO=PRO_
            OCC=OCC_
            KPOINTS=KPOINTS_
        else:
            print(KPOINTS.shape,KPOINTS_.shape)
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

if setfermi==True:
    #EIG_f=EIG
    #EIG_f[OCC<0.01]=-np.inf
    #Nb=np.prod(EIG_f.shape[1:])
    fermi=np.max(EIG[OCC>0.001])
    print(fermi)
    #fermi_ax=np.argmax(EIG_.reshape((2,Nb)))
    print("Fermi is %f" % fermi)
    EIG=EIG-fermi

#EIG=EIG-4.5043

L_INX=[(not i=="") for i in symbol]
if kpath.shape[0]==1:
    PRO=np.concatenate((PRO,PRO),axis=1)
    EIG=np.concatenate((EIG,EIG),axis=1)
    kpath=np.array([0,1])
    intx=np.linspace(0,np.max(kpath),101)
    fig=plt.figure(figsize=(3*ISPIN,10))
else:
    intx=np.linspace(0,np.max(kpath),1001)
    fig=plt.figure(figsize=(8*ISPIN,10))

Elim=(-2-0.17,5-0.17)
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
        #ax_cbar = divider.append_axes('right', size="10%", pad=0.02)
        sc=ax.scatter(np.array([intx.T]*(EIG.shape[2])).T,E,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='Greys')
        for ele in [["S",'y','s'],["Cu",'r','>'],["Sn",'g','v'],["Zn",'b','^']]:
            L_INX=[(i==ele[0]) for i in symbol]
            L1=np.sum(PRO[ispin,:,:,L_INX],axis=0)/np.sum(PRO[ispin],axis=2)
            L=L1#/(L1+L2)a
            intx2=np.linspace(0,np.max(kpath),kpath.shape[0]*2)
            c=griddata(kpath,np.array(L),intx2)
            E=griddata(kpath,EIG[ispin,:],intx2)
            ax.scatter(np.array([intx2.T]*(EIG.shape[2])).T,E,color='',edgecolor=ele[1],s=c*50,marker=ele[2],label=ele[0])
        ax.legend(loc='upper right')
        if any(kpt_label!=""):
            ax.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
            ax.set_xticks(kpath[kpt_label!=""])
            ax.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
        else:
            ax.set_xticks([])
        ax.set_xlim([np.min(kpath),np.max(kpath)])
        ax.set_ylim([Elim[0],Elim[1]])
        ax.set_yticklabels(ax.get_yticks(),fontsize=20)
        #cb=plt.colorbar(sc,cax=ax_cbar,orientation='vertical')
        #cb.set_ticks([0,1])
        #cb.ax.tick_params(labelsize=20)
        #cb.set_ticklabels(['','S'])

    else:
        gs=GridSpec(1,9)
        if ispin==0:
            ax1=plt.subplot(gs[0,0:4])
            sc1=ax1.scatter(np.array([intx.T]*(EIG.shape[2])).T,E,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='jet')
            if any(kpt_label!=""):
                ax1.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
                ax1.set_xticks(kpath[kpt_label!=""])
                ax1.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
            else:
                ax1.set_xticks([])
            ax1.set_xlim([np.min(kpath),np.max(kpath)])
            ax1.set_ylim([Elim[0],Elim[1]])
            ax1.set_yticklabels(ax1.get_yticks(),fontsize=20)
            print("0")
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
            ax2.set_ylim([Elim[0],Elim[1]])
            ax_cbar=plt.subplot(gs[:,8])
            cb=plt.colorbar(sc2,cax=ax_cbar,orientation='vertical')
            cb.set_ticks([0,1])
            cb.ax.tick_params(labelsize=20)
            cb.set_ticklabels(['PTI','Pt'])
        plt.subplots_adjust(wspace=0.02)
plt.tight_layout()
plt.savefig('band.png',dpi=360)
plt.show()

