#!/usr/bin/python3
import ase.io
import numpy as np
import pandas as pd
import re
import matplotlib
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
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
class plot_band:
    def __init__(self):
        pass
    #################################
    def from_kpoints(path='.',KPOINTS='KPOINTS',read=False,save=True):
        atoms  = ase.io.read(path+'/CONTCAR')
        rep    = atoms.get_reciprocal_cell()
        K_files=glob(path+'/'+KPOINTS+'_[0-9]*')
        #print(K_files)
        if K_files==[]:
            K_files=[path+"/"+KPOINTS]
        K_files.sort()
        #print("reading ",K_files)
        KPATH=make_k_path(K_files,rep)   
        return KPATH
    #################################
    def from_procar(path='.',ISPIN=1,read=False,save=True):
        PROCAR_files=glob(path+'/PROCAR_[0-9]*')
        if PROCAR_files==[]:
            PROCAR_files=[path+"/PROCAR"]
        PROCAR_files.sort()
        project=get_procar(PROCAR_files,ISPIN)
        return project
    #################################
    def get_fermi(EIG,OCC,OCC_E=0.001):
        FERMI=np.max(EIG[OCC>OCC_E])
        #print("Fermi is %f" % fermi)
        return FERMI
    #################################
    def set_fermi(EIG,OCC,OCC_E=0.001,FERMI=None):
        if FERMI is None:
            FERMI=np.max(EIG[OCC>OCC_E])
        #print("Fermi is %f" % FERMI)
        EIG=EIG-FERMI
        #print("set fermi to zero")
        return EIG
    #################################
    def get_symbol(CONTCAR_file='./POSCAR',atoms=None):
        if atoms is None:
            try:
                atoms  = ase.io.read(CONTCAR_file)
            except:
                raise("need POSCAR")
        SYMBOL=np.array(atoms.get_chemical_symbols())
        return SYMBOL
    
    def plot_band( kpoints,
                   eigcar,
                   ax=None,
                   Elim=None,
                   ispin=0
                   ):
        if ax is None:
            ax=plt.subplot()
        iE=eigcar.merge(kpoints,how='outer')
        sE=iE[['kpath','eig','band']].drop_duplicates( keep='first', inplace=False)
        sE=sE.sort_index(by=['band','kpath'])
        for i in sE['band'].drop_duplicates( keep='first', inplace=False):
            line=mlines.Line2D(sE[sE['band']==i]['kpath'],sE[sE['band']==i]['eig'],color='black')
            line.set_zorder(0)
            ax.add_line(line)
        return ax
    
    def plot_sp_kline( kpath,
                   ax=None,
                   ):
        if ax is None:
            ax=plt.subplot()
        #print(kpath,[ (i != "")&( not "!" in i) for i in kpath[1]])
        KPATH_SP=kpath[0][[ (i != "")&( not "!" in i) for i in kpath[1]]]
        if len(KPATH_SP) > 0:
            ax.plot([KPATH_SP]*2,np.array([-1000,1000]),color='black')
            ax.set_xticks(KPATH_SP)
            ax.set_xticklabels(kpath[1][[ (i != "")&( not "!" in i) for i in kpath[1]]])#)#,fontsize=24)
       # if sE['kpath'].min()==kpoints['kpath'].max():
        #    ax.set_xlim([-0.5,0.5])
         #   ax.set_xticklabels([])
          #  ax.set_xticks([])
        ##
        return ax

    
    def plot_pband( 
                    plot_type,
                    kpath,
                    eig,
                    pro_group,
                    group_info,
                    ax=None,
                    ispin=0,
                    fermi=0,
                    ):
        ####
        #plt.rc('font',family='Times New Roman')
        #matplotlib.rcParams['mathtext.fontset'] = 'stix'
        #
        divider=make_axes_locatable(plt.gca())
        eig=eig-fermi
        if ax is None:
            ax=plt.subplot()
            #if plot_type!=2:
            # set special K-points
        #ax.set_ylim((Elim[0],Elim[1]))
        #if kpoints['kpath'].min()==kpoints['kpath'].max():
        #    ax.set_xlim([-0.5,0.5])
        #    ax.set_xticklabels([])
        #else:
        #    ax.set_xlim([kpoints['kpath'].min(),kpoints['kpath'].max()])
        #ax.set_yticklabels(ax.get_yticks())#)#,fontsize=24)
        ## plot PBAND  ##
        ##
        #
        # plot use a colorbar
        if plot_type==2:
            if kpath[0].max()==kpath[0].min():
                x=np.linspace(-0.5,0.5,100)
                y=griddata(kpath[0],eig,x)
                c=griddata(kpath[0],pro_group[0],x)
                ax_cbar = divider.append_axes('right', size="50%", pad=0.1)
            else:
                x=np.linspace(kpath[0].min(),kpath[0].max(),1000)
                y=griddata(kpath[0],eig,x)
                c=griddata(kpath[0],pro_group[0],x)
#            sc=plt.scatter(x=kpath[0],y=kp['eig'],c=c,vmin=0,vmax=1,s=2,linewidths=None,marker='o',cmap='jet')
                ax_cbar = divider.append_axes('right', size="5%", pad=0.02)
            #print(x.shape,y.shape,c.shape)
            sc=ax.scatter(x=x.repeat(y.shape[1]),y=y,c=c,vmin=0,vmax=1,linewidths=None,marker='o',cmap='jet')
            cb=plt.colorbar(sc,cax=ax_cbar,orientation='vertical')
            cb.set_ticks([0,1])
            cb.ax.tick_params(labelsize=24)
            print(group_info)
            cb.set_ticklabels([group_info[0],group_info[1]])
     #   plot use symbols
        if plot_type==1:
            #EIG1=EIG[ispin]-1/10
            LS=0
            Tex=[]
            T=[]
            EIG_up1=eig.copy()
            EIG_down1=eig.copy()     
            print(pro_group.shape)
            for i,igroup in enumerate(pro_group):
                label=r''+group_info[i][0]
                EIG_up2=EIG_up1+pro_group[i]/pro_group.max()/20
                EIG_down2=EIG_down1-pro_group[i]/pro_group.max()/20
                print(i)
                if kpath[0].max()==kpath[0].min():
                    for iband in range(eig.shape[-1]):
                        axfill=plt.fill_between([-0.5,0.5],EIG_up1[:,iband].repeat(2),EIG_up2[:,iband].repeat(2),color=group_info[i][1],alpha=0.75,label=label,lw=0)
                        axfill=plt.fill_between([-0.5,0.5],EIG_down1[:,iband].repeat(2),EIG_down2[:,iband].repeat(2),color=group_info[i][1],alpha=0.75,label=label,lw=0)

                        label=None
                else:
                    for iband in range(eig.shape[-1]):
                        axfill=plt.fill_between(kpath[0],EIG_up1[:,iband],EIG_up2[:,iband],color=group_info[i][1],alpha=0.75,label=label,lw=0)
                        axfill=plt.fill_between(kpath[0],EIG_down1[:,iband],EIG_down2[:,iband],color=group_info[i][1],alpha=0.75,lw=0)
                        label=None
                EIG_up1=EIG_up2
                EIG_down1=EIG_down2
            ax.legend(loc='upper right')#)#,fontsize=24)
            ax_cbar=None
            divider=make_axes_locatable(plt.gca())
        return ax,ax_cbar,divider


