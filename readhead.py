#!/usr/bin/python3
import sys
from glob import glob

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
from readvasp import *
import plotting
iniset={}
class plotset(object):
    def __init__(self):
        pass
    def group_info_from_input(self,iniset,type):
        group_tag=[]
        group_info=[]
        if iniset['plot_type']%10==0:
            pass
        elif (type==2) | (iniset['plot_type']%10==1):
            tmp=(input("Element or atomic NO., specified orbital with \"_\" ,split with \",\", label color ;\n e.g : Cu,Ag_py,1-3,6_py,8-10_d CuAnTag red\n empty for end\n"))
            while not tmp == "":
                type1=tmp.split()
                group_tag.append(type1[0])
                group_info.append(type1[1:])
                tmp=(input())
        elif (type==1) & (iniset['plot_type']%10==2):
            tmp=(input("Element or atomic NO., specified orbital with \"_\" ,split with \",\", label;\n e.g : Cu,Ag_py,1-3,6_py,8-10_d CuAnTag\n empty for end\n partA:\n"))
            group_tag.append(tmp.split()[0])
            group_info.append(tmp.split()[1])
            tmp=input("partB:\n")
            if not tmp =="":
                type1=tmp.split()
                group_tag.append(type1[0])
                group_info.append(type1[1])
            else:
                group_info.append("")
        return group_tag,group_info
    def Elim_from_input(self):
        tmp=input("Energy range for y axis,split with \",\"\n e.g -2,5  default is [-2,5]\n")
        if tmp=="":
            Elim=(-2,5)
        else:
            tmp=tmp.split()
            Elim=(float(tmp[0]),float(tmp[1]))
        return Elim
    def tdm_from_input(self):
        if input("Wether to plot tdm of band. (T/F) \n").lower() == 't':
            tmp=input("tdm band, e.g 1,6 1-2,6-9  \n")
            tdm=np.empty(shape=[0,2])
            while not tmp =="":
                tttdm1,tttdm2=tmp.split()
                if '-' in tttdm1:
                    ttdm1=np.arange(int(tttdm1.split('-')[0])-1,int(tttdm1.split('-')[1]))
                else:
                    ttdm1=np.array(tttdm1,dtype=np.int)-1
                if '-' in tttdm2:
                    ttdm2=np.arange(int(tttdm2.split('-')[0])-1,int(tttdm2.split('-')[1]))
                else:
                    ttdm2=np.array(tttdm2,dtype=np.int)-1
                tdm1,tdm2=np.meshgrid(ttdm1,ttdm2)
                ttdmlist=np.array([tdm1.flatten(),tdm2.flatten()]).reshape(2,-1).T
                tdm=np.append(tdm,ttdmlist).reshape(-1,2)
                tmp=input("tdm band \n")
            tdmlist_s=tdm.tolist()
            tdmlist_s.sort()
            it = itertools.groupby(tdmlist_s)
            tdmlist=[]
            for k, g in it:
                tdmlist.append(k)
            else:
                tdm=False
            tdmlist=np.array(tdmlist)
            tdm=tdmlist[tdmlist[:,0]<tdmlist[:,1]]
        else:
            tdm=False
        return tdm
    def figset_from_input(self,type):
        tmp=input("Figsize, split with \",\"\n e.g 6,8 default is [6,8]\n")
        if tmp =="":
            if type==1:
                figsize=(6,8)
            elif type==2:
                figsize=(8,4)
        else:
            tmp=tmp.split()
            figsize=(float(tmp[0]),float(tmp[1]))
        tmp=input("fontsize, e.g 24 default is 24 \n")
        if tmp =="":
            fontsize=24
        else:
            fontsize=float(tmp)
        tmp=input("Fonttype, T: Times New Roman, F: Default \n")
        plt.rc('font',family='serif',size=fontsize)
        if tmp.lower()=='T':
            plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
            matplotlib.rcParams['mathtext.fontset'] = 'stix'
        return figsize,fontsize
    def L_fermi_from_input(self):
        tmp=input(  "Where to set as fermi? (T/VBM,CBM,band,Energy/F/0)\n"
                    "F/0 will not change anything.\n" 
                    "e.g VBM/CBM/band_18/-2.5/F/0\n"
                    "default is True")
        if (tmp.lower() == 't') | (tmp==""):
            iniset['L_fermi'] = True
            iniset['specific_fermi'] = False
        elif (tmp.lower() == 'f') :
            iniset['L_fermi'] = False
        else:
            iniset['L_fermi'] = True
            iniset['specific_fermi'] = True
            iniset['specific_fermi_value'] = float(tmp)
        return iniset
    def read_from_input(self,type):
        if type==1:
            #iniset['I_spin']=int(input( "spin, e.g 1 or 2 \n"))
            iniset['plot_type']=int(input(  "plot type, \n"
                                            "e.g : nomal_band/0,fat_band/1,colormap_band/2\n"))
        elif type==2:
            iniset['plot_type']=int(input(  "plot type, \n"
                                            "e.g : total dos/0,pdos/2\n"))
        iniset['group_tag'],iniset['group_info']=self.group_info_from_input(iniset,type)
        iniset['Elim']=self.Elim_from_input()
        if type==1:
            iniset['tdm']=self.tdm_from_input()
        iniset['L_fermi']=self.L_fermi_from_input()
        iniset['figsize'],iniset['fontsize']=self.figset_from_input(type)
        iniset['intr']=1000
        return iniset

def plotband():
    iniset=plotset().read_from_input(1)
    if iniset['plot_type']==0:
        data=from_eigenval.get_eigenvalue()
    else:
        data=from_procar.get_procar()
    kpt=from_kpoints.get_kpoints()
    poscar=from_poscar.get_poscar('./POSCAR')
    symbollist=poscar.get_symbollist()
    kpt.make_k_path(kpt.kpoint,kpt.sp_kpt_label,poscar.bcell)
    if (iniset['plot_type']==0):
        if (data.I_spin>1):
            tmp=input('in one figure? T/F\n')
            if 't' in tmp.lower():
                iniset['onefigure']=True
            else:
                iniset['onefigure']=False
    else:
        iniset['onefigure']=False
    for ispin in range(data.I_spin):
        if (not iniset['onefigure'])  | (ispin==0):
            fig=plt.figure(figsize=iniset['figsize'])
            ax=plt.subplot()
        ax.set_ylim(iniset['Elim'])
        pltband=plotting.plot_band(ispin=ispin,kpath=kpt.kpath,sp_kpt_label=kpt.sp_kpt_label,data=data,iniset=iniset)
        if iniset['plot_type'] > 0:
            ax,divider,ax_cbar=pltband.plot_band(   plot_type=iniset['plot_type'],
                                                    pro_group=data.set_group(iniset['group_tag'],symbollist)[:,ispin],
                                                    group_info=iniset['group_info'],
                                                    ax=ax)
        else:
            ax,divider,ax_cbar=pltband.plot_band(   iniset['plot_type'],ax=ax)
        ax.set_xlim((pltband.kpath_intr.min(),pltband.kpath_intr.max()))
        fig.subplots_adjust(0.3,0.2,0.5,0.9)
        ax.legend(loc='right', bbox_to_anchor=(3.5, 0.6),frameon=False)
        if iniset['tdm']!=False:
            tdmeig=np.zeros((data.I_spin,kpt.kpath.size,iniset['tdm'].size))
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax_tdm= divider.append_axes('bottom', size="15%", pad=0.1)
            wfc=from_wavecar.get_wavecar()
            tdm=wfc.get_tdm()
            tdm.fermi=0
            for i,itdm in enumerate(iniset['tdm']):
                tdmeig[ispin,:,i]=tdm.calc_tdm_in_band(np.arange(len(kpt.kpath),dtype=int)+1,ispin+1,int(itdm[0])+1,ispin+1,int(itdm[1])+1,norm=True)[:,-1].real
            tdm.eig=tdmeig
            plttmdband=plotting.plot_band(ispin=ispin,kpath=kpt.kpath,sp_kpt_label=kpt.sp_kpt_label,data=tdm,iniset=iniset)
            ax_tdm,_,_=plttmdband.plot_band(plot_type=0,ax=ax_tdm,divider=divider)
            ax_tdm.set_xlim((pltband.kpath_intr.min(),pltband.kpath_intr.max()))
            ax_tdm.set_ylabel(r'$\mathrm{TDM\ [ Debyte^2 ]}$')
            #ax_tdm.set_yticklabels(ax.get_yticks())#,fontsize=24)
            ax_tdm.set_ylim((0,np.max(tdm.eig)))
            #plt.legend( markerscale=0.5)#,fontsize=24)
        #plt.tight_layout()
        #plt.savefig('band_%d.png' % ispin,dpi=300)
    plt.show()

def plotdos():
    iniset=plotset().read_from_input(2)
    data=from_doscar.get_doscar()
    plt.figure(figsize=iniset['figsize'])
    ax=plt.subplot()
    ax.set_xlim(iniset['Elim'])
    pltdos=plotting.plot_dos(data,iniset=iniset)
    if iniset['plot_type']>0:
        poscar=from_poscar.get_poscar('./POSCAR')
        symbollist=poscar.get_symbollist()
        ax=pltdos.plot_dos(iniset['plot_type'],data.set_group(iniset['group_tag'],symbollist),iniset['group_info'],ax=ax)
    else:
        ax=pltdos.plot_dos(iniset['plot_type'],ax=ax)
    Ymax=data.total[:, ( data.Energy>iniset['Elim'][0] ) & ( data.Energy<iniset['Elim'][1] ) ].max()*1.2
    if data.I_spin==2:
        Ymin=-Ymax
    else:
        Ymin=0
    ax.set_ylim((   Ymin,   Ymax))
    plt.show()

if __name__ == "__main__":
    tmp=int(input('Band:[1]/DOS[2]/CHARGE[3]\n'))
    if tmp==1:
        plotband()
    elif tmp==2:
        plotdos()
