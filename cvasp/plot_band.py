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
import seaborn as sns
class plot_band:
    def __init__(self):
        pass
    #################################
    def from_kpoints(path='.',KPOINTS='KPOINTS',read=False,save=True):
        atoms  = ase.io.read(path+'/CONTCAR')
        rep    = atoms.get_reciprocal_cell()
        K_files=glob(path+'/'+KPOINTS+'_[0-9]*')
        print(K_files)
        if K_files==[]:
            K_files=[path+"/"+KPOINTS]
        K_files.sort()
        #print("reading ",K_files)
        kpoints=make_k_path(K_files,rep)   
        return kpoints

        
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
    def add_info_symbol(project,CONTCAR_file='./CONTCAR',atoms=None):
        if atoms is None:
            try:
                atoms  = ase.io.read(CONTCAR_file)
            except:
                raise("need CONTCAR")
        symbol=pd.DataFrame( )
        symbol['ion']= np.arange(atoms.get_number_of_atoms())+1
        symbol['symbol']= np.array(atoms.get_chemical_symbols())
        project=project.merge(symbol,how='outer')
        return project
    
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
    
    def plot_sp_kline( kpoints,
                   ax=None,
                   ):
        if ax is None:
            ax=plt.subplot()
        iE=kpoints
        sE=iE[(iE['kpt_label']!="")&(iE['kpt_label']!="!")][['kpath','kpt_label']].drop_duplicates( keep='first', inplace=False)
        if len(sE) > 0:
            ax.plot([sE['kpath']]*2,np.array([-1000,1000]),color='black')
            ax.set_xticks(sE['kpath'])
            ax.set_xticklabels(sE['kpt_label'],fontsize=24)
        if sE['kpath'].min()==kpoints['kpath'].max():
            ax.set_xlim([-0.5,0.5])
            ax.set_xticklabels([])
            ax.set_xticks([])
        ##
        return ax

    
    def plot_pband( 
                    plot_type,
                    kpoints,
                    project,
                    ax=None,
                    Elim=None,
                    ispin=0,
                    type2_Ele_A="",
                    type2_Ele_B="",
                    type1_Ele=False,
                    set_fermi=True
                    ):
        ####
        if set_fermi==True:
            project['eig']=project['eig']-project[project['occ']>0.001]['eig'].max()
        vbm=project[project['occ']>0.001]['eig'].max()
        cbm=project[project['occ']<=0.001]['eig'].min()
        promap={'s':['s'],'py':['py'],'pz':['pz'],'px':['px'],'dxy':['dxy'],'dyz':['dyz'],'dz2':['dz2'],'dx2-y2':['dx2-y2'],'dxz':['dxz'],'p':['px','py','pz'],'d':['dxy','dyz','dz2','dx2-y2','dxz']}
        protex={'s':'s','py':'p_y','pz':'p_z','px':'p_x','dxy':'d_{xy}','dyz':'d_{yz}','dz2':'d_{z^2}','dx2-y2':'d_{x^2-y^2}','dxz':'d_{xz}','p':'p','d':'d'}
        #plt.rc('font',family='Times New Roman')
        #matplotlib.rcParams['mathtext.fontset'] = 'stix'
        iE=project[project['spin']==ispin]
        minband=iE[(iE['eig']>=Elim[0])&(iE['eig']<=Elim[1])]['band'].min()
        maxband=iE[(iE['eig']>=Elim[0])&(iE['eig']<=Elim[1])]['band'].max()
        iE=iE[(iE['band']>=minband)&(iE['band']<=maxband)]
        iE=iE.merge(kpoints,how='outer')
        #
        if ax is None:
            ax=plt.subplot()
            #if plot_type!=2:
            # set special K-points
        ax.set_ylim((Elim[0],Elim[1]))
        if kpoints['kpath'].min()==kpoints['kpath'].max():
            ax.set_xlim([-0.5,0.5])
            ax.set_xticklabels([])
        else:
            ax.set_xlim([kpoints['kpath'].min(),kpoints['kpath'].max()])
        ax.set_yticklabels(ax.get_yticks(),fontsize=24)
        ## plot PBAND  ##
        ##
        #
        tE=iE[['tot','band','kpt']].groupby(['kpt','band']).sum().reset_index()
        tE.rename(columns={'tot':'Ttot'}, inplace = True)
        iE=iE.merge(tE,how='outer')
        kp=iE[['kpath','kpt','band','eig','Ttot']].drop_duplicates()
        # plot use a colorbar
        if plot_type==2:
            A=[]
            TexA=[] 
            for iA in type2_Ele_A:
                iAs=iA.split()[0].split('_')
                Ele_A=iAs[0]
                if len(iAs)>1:
                    orbit=promap[iAs[1]]
                    TexA.append(r'$\mathrm{'+Ele_A+'_{'+protex[iAs[1]]+'}}$')
                else:
                    TexA.append(r'$\mathrm{'+Ele_A+'}$')
                    orbit=['tot']
                mE=iE[(iE['symbol']==Ele_A)][orbit+['band','kpt']].groupby(['kpt','band']).sum().reset_index()
                mE[iA]=mE[orbit].apply(lambda x: x.sum(), axis=1)
                kp=kp.merge(mE[['band','kpt',iA]],how='outer')
                A.append(iA)
               #L_INX_A=[i in type2_Ele_A for i in symbol]
            kp['totA']=kp[A].apply(lambda x: x.sum(), axis=1)
            if type2_Ele_B:
                B=[]
                TexB=[]
                for iB in type2_Ele_B:
                    iBs=iB.split()[0].split('_')
                    Ele_B=iBs[0]
                    if len(iBs)>1:
                        orbit=promap[iBs[1]]
                        TexB.append(r'$\mathrm{'+Ele_B+'_{'+protex[iBs[1]]+'}}$')
                    else:
                        orbit='tot'
                    mE=iE[(iE['symbol']==Ele_B)][orbit+['band','kpt']].groupby(['kpt','band']).sum().reset_index()
                    mE[iB]=mE[orbit].apply(lambda x: x.sum(), axis=1)
                    kp=kp.merge(mE[['band','kpt',iB]],how='outer')
                    B.append(iB)
                   #L_INX_A=[i in type2_Ele_A for i in symbol]
                kp['totB']=kp[B].apply(lambda x: x.sum(), axis=1)
                c=kp['totA']/(kp['totB']+kp['totA'])
            else:
                c=kp['totA']/kp['Ttot']
            divider=make_axes_locatable(plt.gca())
            if kp['kpath'].max()==kp['kpath'].min():
                xx=np.linspace(-0.5,0.5,100)
                x=np.array([])
                y=np.array([])
                c=np.array([])
                for ikp in kp.groupby(['band']):
                    x=np.append(x,xx)
                    y=np.append(y,ikp[1]['eig'].values.repeat(100))
                    c=np.append(c,ikp[1]['totA'].values/ikp[1]['Ttot'].values.repeat(100))
                ax_cbar = divider.append_axes('right', size="50%", pad=0.1)
            else:
                xx=np.linspace(kp['kpath'].min(),kp['kpath'].max(),1000)
                x=np.array([])
                y=np.array([])
                c=np.array([])
                for ikp in kp.groupby(['band']):
                    x=np.append(x,xx)
                    y=np.append(y,griddata(ikp[1]['kpath'].values,ikp[1]['eig'].values,xx))
                    c=np.append(c,griddata(ikp[1]['kpath'].values,ikp[1]['totA'].values/ikp[1]['Ttot'].values,xx))
#            sc=plt.scatter(x=kp['kpath'],y=kp['eig'],c=c,vmin=0,vmax=1,s=2,linewidths=None,marker='o',cmap='jet')
                ax_cbar = divider.append_axes('right', size="5%", pad=0.02)
            sc=ax.scatter(x=x,y=y,c=c,vmin=0,vmax=1,s=2,linewidths=None,marker='o',cmap='jet')
            cb=plt.colorbar(sc,cax=ax_cbar,orientation='vertical')
            cb.set_ticks([0,1])
            cb.ax.tick_params(labelsize=24)
            cb.set_ticklabels(['',",".join(TexA)])
     #   plot use symbols
        if plot_type==1:
            #EIG1=EIG[ispin]-1/10
            LS=0
            Tex=[]
            T=[]
            kp['T11']=kp['eig']
            kp['T21']=kp['eig']
            for ele in type1_Ele:
                TexA=[]
                A=[]
                for iA in ele[0]:
                    iAs=iA.split()[0].split(r'_')
                    Ele_A=iAs[0]
                    if len(iAs)>1:
                        orbit=promap[iAs[1]]
                        TexA.append(r'$\mathrm{'+Ele_A+'_{'+protex[iAs[1]]+'}}$')
                    else:
                        TexA.append(r'$\mathrm{'+Ele_A+'}$')
                        orbit=['tot']
                    mE=iE[(iE['symbol']==Ele_A)][orbit+['band','kpt']].groupby(['kpt','band']).sum().reset_index()
                    mE[iA]=mE[orbit].apply(lambda x: x.sum(), axis=1)
                    kp=kp.merge(mE[['band','kpt',iA]],how='outer')
                    A.append(iA)
                Tex.append(",".join(TexA))
                label=",".join(ele[0])
                kp[label]=kp[A].apply(lambda x: x.sum(), axis=1)
                width=100/(Elim[1]-Elim[0])
                kp['T12']=kp['T11']+kp[label]/width
                kp['T22']=kp['T21']-kp[label]/width
                if kp['kpath'].max()==kp['kpath'].min():
                    for i in range(kp['band'].min(),kp['band'].max()+1):
                        kkp=kp[kp['band']==i]
                        axfill=plt.fill_between([-0.5,0.5],kkp['T11'].values.repeat(2),kkp['T12'].values.repeat(2),color=ele[1],alpha=0.75,label=label)
                        axfill=plt.fill_between([-0.5,0.5],kkp['T21'].values.repeat(2),kkp['T22'].values.repeat(2),color=ele[1],alpha=0.75)
                        label=None
                else:
                    for i in range(kp['band'].min(),kp['band'].max()+1):
                        kkp=kp[kp['band']==i]
                        axfill=plt.fill_between(kkp['kpath'],kkp['T11'],kkp['T12'],color=ele[1],alpha=0.75,label=label)
                        axfill=plt.fill_between(kkp['kpath'],kkp['T21'],kkp['T22'],color=ele[1],alpha=0.75)
                        label=None
                kp['T11']=kp['T12']
                kp['T21']=kp['T22']
            ax.legend(loc='upper right',fontsize=24)
            ax_cbar=None
            divider=None

        if plot_type==1:
            #EIG1=EIG[ispin]-1/10
            LS=0
            Tex=[]
            T=[]
            kp['T11']=kp['eig']
            kp['T21']=kp['eig']
            for ele in type1_Ele:
                TexA=[]
                A=[]
                for iA in ele[0]:
                    iAs=iA.split()[0].split(r'_')
                    Ele_A=iAs[0]
                    if len(iAs)>1:
                        orbit=promap[iAs[1]]
                        TexA.append(r'$\mathrm{'+Ele_A+'_{'+protex[iAs[1]]+'}}$')
                    else:
                        TexA.append(r'$\mathrm{'+Ele_A+'}$')
                        orbit=['tot']
                    mE=iE[(iE['ion']==Ele_A)][orbit+['band','kpt']].groupby(['kpt','band']).sum().reset_index()
                    mE[iA]=mE[orbit].apply(lambda x: x.sum(), axis=1)
                    kp=kp.merge(mE[['band','kpt',iA]],how='outer')
                    A.append(iA)
                Tex.append(",".join(TexA))
                label=",".join(ele[0])
                kp[label]=kp[A].apply(lambda x: x.sum(), axis=1)
                width=100/(Elim[1]-Elim[0])
                kp['T12']=kp['T11']+kp[label]/width
                kp['T22']=kp['T21']-kp[label]/width
                if kp['kpath'].max()==kp['kpath'].min():
                    for i in range(kp['band'].min(),kp['band'].max()+1):
                        kkp=kp[kp['band']==i]
                        axfill=plt.fill_between([-0.5,0.5],kkp['T11'].values.repeat(2),kkp['T12'].values.repeat(2),color=ele[1],alpha=0.75,label=label)
                        axfill=plt.fill_between([-0.5,0.5],kkp['T21'].values.repeat(2),kkp['T22'].values.repeat(2),color=ele[1],alpha=0.75)
                        label=None
                else:
                    for i in range(kp['band'].min(),kp['band'].max()+1):
                        kkp=kp[kp['band']==i]
                        axfill=plt.fill_between(kkp['kpath'],kkp['T11'],kkp['T12'],color=ele[1],alpha=0.75,label=label)
                        axfill=plt.fill_between(kkp['kpath'],kkp['T21'],kkp['T22'],color=ele[1],alpha=0.75)
                        label=None
                kp['T11']=kp['T12']
                kp['T21']=kp['T22']
            ax.legend(loc='upper right',fontsize=24)
            ax_cbar=None
            divider=None

        return ax,ax_cbar,divider


