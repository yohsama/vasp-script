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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines
class plot_band:
    def __init__(self):
        pass
    #################################
    def from_kpoints(path='.',read=False,save=True):
        atoms  = ase.io.read(path+'/CONTCAR')
        rep    = atoms.get_reciprocal_cell()
        if read:
            try:
                kpath     = np.load(path+'/kpath.npy')
                kpt_label = np.load(path+'/kpt_label.npy')
            except:
                K_files=glob(path+'/KPOINTS_*')
                if K_files==[]:
                    K_files=[path+"/KPOINTS"]
                K_files.sort()
                #print("reading ",K_files)
                kpath,kpt_label,KPOINTS=make_k_path(K_files,rep)  
        else:
            K_files=glob(path+'/KPOINTS_*')
            if K_files==[]:
                K_files=[path+"/KPOINTS"]
            K_files.sort()
            #print("reading ",K_files)
            kpath,kpt_label,KPOINTS=make_k_path(K_files,rep)   
        if save:
            np.save(path+'/kpath.npy',kpath)
            np.save(path+'/kpt_label.npy',kpt_label)
        return kpath,kpt_label,KPOINTS

        
    #################################
    def from_procar(path='.',ISPIN=1,read=False,save=True):
        if read:
            try:
                PRO       = np.load(path+'/PRO.npy')
                EIG       = np.load(path+'/EIG.npy')
                OCC       = np.load(path+'/OCC.npy')
                KPOINTS   = np.load(path+'/KPOINTS.npy')
            except:
                print("read fail")
                i=0
                PROCAR_files=glob(path+'/PROCAR_*')
                if PROCAR_files == []:
                    PROCAR_files=[path+"/PROCAR"]
                EIG,OCC,PRO,KPOINTS=get_procar(PROCAR_files,ISPIN)
        else:
            i=0
            PROCAR_files=glob(path+'/PROCAR_*')
            if PROCAR_files == []:
                PROCAR_files=[path+"/PROCAR"]
            EIG,OCC,PRO,KPOINTS=get_procar(PROCAR_files,ISPIN)
        if save:    # save dat for read
            np.save(path+'/PRO.npy',PRO)
            np.save(path+'/EIG.npy',EIG)
            np.save(path+'/OCC.npy',OCC)
            np.save(path+'/KPOINTS.npy',KPOINTS)
        return EIG,OCC,PRO,KPOINTS
    #################################
    def get_fermi(EIG,OCC,OCC_E=0.001):
        fermi=np.max(EIG[OCC>OCC_E])
        print("Fermi is %f" % fermi)
        return fermi
    #################################
    def set_fermi(EIG,OCC,OCC_E=0.001):
        fermi=np.max(EIG[OCC>OCC_E])
        print("Fermi is %f" % fermi)
        EIG=EIG-fermi
        print("set fermi to zero")
        return EIG
    #################################
    def plot_band(kpath,
                    kpt_label,
                    EIG,
                    Elim=(-2,5),
                    PRO=False,
                    figsize=False,
                    intd=False,
                    atoms=False,
                    type2_Ele_A="",
                    type2_Ele_B="",
                    type1_Ele=False,
                    type2_Number_A=False,
                    type2_Number_B=False,
                    plot_type=1
                    ):
        ####
        print(kpath.shape)
        if kpath.shape[0]==1: # Gamma only
            PRO=np.concatenate((PRO,PRO),axis=1)
            EIG=np.concatenate((EIG,EIG),axis=1)
            kpath=np.array([0,0.1])
        ####
        #
        #
        # reduce the size of EIG/PRO by Elim
        EIG_=np.max(np.max(EIG,axis=0),axis=0)
        RL=np.argmax(EIG_>Elim[0])
        EIG_=np.min(np.min(EIG,axis=0),axis=0)
        RH=np.argmax(EIG_>Elim[1])
        if RH==0:
            RH=EIG.shape[2]
        PRO=PRO[:,:,RL:RH]
        EIG=EIG[:,:,RL:RH]
        print(RH,RL)
        ##
        #
        ISPIN=EIG.shape[0]
        ##   set size of image
        if not figsize:
            if kpath.shape[0]>2:
                width=6
            else:
                width=2
            hight=8
            figsize=(((ISPIN)*width+0.5*np.max((len(type2_Ele_A),len(type2_Ele_B))),hight))
        else:
            width=(int(figsize[0]-0.5*np.max((len(type2_Ele_A),len(type2_Ele_B))))/2/ISPIN)
        fig=plt.figure(figsize=figsize)
        ##
        #
        if not atoms:
            atoms  = ase.io.read('./CONTCAR')
            pos    = atoms.get_scaled_positions()
            symbol = atoms.get_chemical_symbols()
            rep    = atoms.get_reciprocal_cell() 
        if not intd:
            intd=np.max(kpath)/0.01
        x=np.linspace(0,np.max(kpath),intd*10)
        if ISPIN==2:
            gs=GridSpec(1, 3, width_ratios=[3, 3,1])
        else:
            gs = GridSpec(1, 2, width_ratios=[3, 1])
        for ispin in range(ISPIN):
            print(EIG[ispin].shape)
            y=griddata(kpath,EIG[ispin,:],x) 
            print(int(ispin*width*2),int(2*width*(ispin+1)))
            ax=plt.subplot(gs[0,ispin])
            #ax=plt.subplot(gs[0,int(ispin*width*2):int((ispin+1)*width*2)])
            ####   plot regular band image
            #
            for i_y in y.T:
                line=mlines.Line2D(x,i_y,color='black')
                line.set_zorder(0)
                ax.add_line(line)
            # set special K-points
            if any(kpt_label!=""):
                    ax.plot([kpath[kpt_label!=""].T]*2,np.array([-100,100]),color='black')
                    ax.set_xticks(kpath[kpt_label!=""])
                    ax.set_xticklabels(kpt_label[kpt_label!=""],fontsize=20)
            else:
                    ax.set_xticks([])
            ##
            ax.set_xlim([np.min(kpath),np.max(kpath)])
            ax.set_ylim([Elim[0],Elim[1]])
            ax.set_yticklabels(ax.get_yticks(),fontsize=20)
            #
            if ispin==1:
                ax.set_yticks([])
            #
            #########
            ## plot PBAND  ##
            ##
            #
            if PRO.any():
                # plot use a colorbar
                if plot_type==2:
                    if type2_Ele_A:
                        L_INX_A=[i in type2_Ele_A for i in symbol]
                    if type2_Ele_B:
                        L_INX_B=[i in type2_Ele_B for i in symbol]
                    L_A=np.sum(PRO[ispin,:,:,L_INX_A],axis=0)/np.sum(PRO[ispin],axis=2)
                    L=L_A
                    if type2_Ele_B:
                        L_B=np.sum(PRO[ispin,:,:,L_INX_B],axis=0)/np.sum(PRO[ispin],axis=2)
                        L=L_A/(L_A+L_B)
                    c=griddata(kpath,np.array(L),x)
                    sc=ax.scatter(np.array([x.T]*(EIG.shape[2])).T,y,vmin=0, vmax=1,c=c,s=2,linewidths=None,marker='o',cmap='jet')
                    if ispin==ISPIN-1:
                        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                        cbar_ax = inset_axes(ax,
                                   width=0.2, # width = 10% of parent_bbox width
                                   height="100%", # height : 50%
                                   loc=6,
                                   bbox_to_anchor=(1.05, 0., 1, 1),
                                   bbox_transform=ax.transAxes,
                                   borderpad=0,
                               )
                        cb=plt.colorbar(sc,cax=cbar_ax,orientation='vertical')
                        cb.set_ticks([0,1])
                        cb.ax.tick_params(labelsize=20)
                        cb.set_ticklabels(['',",".join(type2_Ele_A)])
                #plot use symbols
                if plot_type==1:
                    for ele in type1_Ele:
                        L_INX=[ i in ele[0] for i in symbol]
                        L=np.sum(PRO[ispin,:,:,L_INX],axis=0)/np.sum(PRO[ispin],axis=2)
                        intx2=np.linspace(0,np.max(kpath),intd)
                        c=griddata(kpath,np.array(L),intx2)
                        E=griddata(kpath,EIG[ispin,:],intx2)
                        print(ele[0])
                        ax.scatter(np.array([intx2.T]*(EIG.shape[2])).T,E,color='',edgecolor=ele[1],s=c*100,marker=ele[2],label=ele[0])
                    if ispin==ISPIN-1:
                        ax.legend(loc='upper right')

#        plt.tight_layout()
        plt.subplots_adjust(left=0.3, right=0.7, top=0.9, bottom=0.1)
#        plt.tight_layout()
        plt.savefig('band.png',dpi=360)
        plt.show()        
        

