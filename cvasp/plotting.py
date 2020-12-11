#!/usr/bin/python3
import os
import re
import sys
from glob import glob
from multiprocessing import pool

import ase.io
import matplotlib
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata


class plot_band(object):
    def __init__(self,ispin,kpath,sp_kpt_label,data,iniset):
        self.kpath=kpath
        self.sp_kpt_label=sp_kpt_label
        self.ispin=ispin
        self.data=data
        if iniset['L_fermi']:
            self.L_fermi=True
            if iniset['specific_fermi']:
                self.fermi=iniset['specific_fermi_value']
            else:
                self.fermi=self.data.fermi
        else:
            self.L_fermi=False
            self.fermi=0
        self._onefigure_=iniset['onefigure']
    
    def plot_band(self,plot_type,intr=1000,pro_group=None,group_info=None,ax=None,divider=None,ax_cbar=None):

        if divider is None:
            divider=make_axes_locatable(plt.gca())
        if ax is None:
            ax=plt.subplot()
        if self.kpath.max()==self.kpath.min():
            """
            if only single kpoints
            expand the data to 100 points width
            """
            self.L_gam=True
            self.kpath_intr=np.linspace(-0.5,0.5,100)
            self.eig_intr=griddata(np.array((-0.5,0.5)),(self.data.eig[self.ispin]-self.fermi).repeat(2,axis=0),self.kpath_intr)
        else:
            self.L_gam=False
            """
            integrate the point, it will makes the type1 - type3 looks like a curve, not a series scattered points.
            """
            self.kpath_intr=np.linspace(self.kpath.min(),self.kpath.max(),intr)
            self.eig_intr=griddata(self.kpath,(self.data.eig[self.ispin]-self.fermi),self.kpath_intr)
           
        """
        type0:
            normal bandgap picture
        type1:
            Project band and marked with various shapes; default is a colored band made up of circles.
        type2:
            colormap of two elements;
        type3:
            colormap of three elements;
        
        You can combine these type by plot multiple times

        """
        if plot_type==0:
            ax,divider,ax_cbar=self.plot_type_0(ax=ax,divider=divider,ax_cbar=ax_cbar)
        elif plot_type==1:
            ax,divider,ax_cbar=self.plot_type_1(pro_group=pro_group,group_info=group_info,ax=ax,divider=divider,ax_cbar=ax_cbar)
        elif plot_type==2:
            ax,divider,ax_cbar=self.plot_type_2(pro_group=pro_group,group_info=group_info,ax=ax,divider=divider,ax_cbar=ax_cbar)
        ax,divider,ax_cbar=self.plot_sp_kline(ax=ax,divider=divider,ax_cbar=ax_cbar)
        if not self.L_fermi:
            ax.set_ylabel(r'$\mathrm{E-E_f (eV)}$')
        else:
            ax.set_ylabel(r'$\mathrm{E (eV)}$')
        return ax,divider,ax_cbar

    def plot_sp_kline(self,ax=None,divider=None,ax_cbar=None):
        if ax is None:
            ax=plt.subplot()
        # "!" is marked for the scf part in HSE/metaGGA calculate and will be exclude in the plot
        self.kpath_sp=self.kpath[[ (i != "")&( not "!" in i) for i in self.sp_kpt_label]] 
        ax.set_xticklabels([])
        if (len(self.kpath_sp) > 0 ) :
            if (len(self.kpath) > 1):
                ax.plot([self.kpath_sp]*2,np.array([-1000,1000]),color='black')
            ax.set_xticks(self.kpath_sp)
            ax.set_xticklabels(self.sp_kpt_label[[ (i != "")&( not "!" in i) for i in self.sp_kpt_label]])
        return ax,divider,ax_cbar


    def plot_type_0(self,ax=None,divider=None,ax_cbar=None):
        if ax is None:
            ax=plt.subplot()  
        if self.data.eig.shape[0] > 1:
            label=['spin up','spin down'][self.ispin]
        else:
            label=None
        if self._onefigure_ & self.ispin == 1:
            c='red'
        else:
            c='black'
        ax.plot(self.kpath_intr,self.eig_intr[:,0],c=c,label=label)
        if self.eig_intr.shape[1]>1:
            ax.plot(self.kpath_intr,self.eig_intr[:,1:],c=c)
        return ax,divider,ax_cbar

    def plot_type_1(self,pro_group,group_info,ax=None,divider=None,ax_cbar=None):
        if ax is None:
            ax=plt.subplot()
        pro_group=np.cumsum(pro_group[::-1],axis=0)[::-1]
        for i,igroup in enumerate(pro_group):
            label=r''+group_info[i][0]
            if self.L_gam:
                prodata=griddata(np.array((-0.5,0.5)),igroup.repeat(2,axis=0),self.kpath_intr)
            else:
                prodata=griddata(self.kpath,igroup,self.kpath_intr)
            ax.scatter(self.kpath_intr.repeat(self.eig_intr.shape[1]),self.eig_intr,c=group_info[i][1],
                        s=prodata*50,linewidths=None,marker='o',cmap='jet',label=label)
        ax.legend(loc='right', bbox_to_anchor=(1, 0.6))
        return ax,divider,ax_cbar

    def plot_type_2(self,pro_group,group_info,ax=None,divider=None,ax_cbar=None):
        if ax is None:
            ax=plt.subplot()
        if pro_group.shape[0] > 1 :
            tsum=pro_group.sum(0)
            tsum[tsum==0]=np.nan
            if self.L_gam:
                tot=griddata(np.array((-0.5,0.5)),pro_group.sum(0).repeat(2,axis=0),self.kpath_intr)
                igroup=griddata(np.array((-0.5,0.5)),(pro_group[0]/tsum).repeat(2,axis=0),self.kpath_intr)
            else:
                tot=griddata(self.kpath,pro_group.sum(0),self.kpath_intr)
                igroup=griddata(self.kpath,pro_group[0]/tsum,self.kpath_intr)
        else:
            if self.L_gam:
                igroup=griddata(np.array((-0.5,0.5)),pro_group[0].repeat(2,axis=0),self.kpath_intr)
            else:
                igroup=griddata(self.kpath,pro_group[0],self.kpath_intr)
            tot=1
        sc=ax.scatter(self.kpath_intr.repeat(self.eig_intr.shape[1]),self.eig_intr,c=igroup,vmin=0,vmax=1,
                        s=tot,linewidths=None,marker='o',cmap='jet')
        cb=plt.colorbar(sc,cax=ax_cbar,orientation='vertical')
        cb.set_ticks([0,1])
        cb.ax.tick_params()
        cb.set_ticklabels([r''+group_info[1],r''+group_info[0]])
        return ax,divider,ax_cbar

class plot_dos(object):
    def __init__(self,dos,iniset):
        self.dos=dos
        if iniset['L_fermi']:
            self.L_fermi=True
            if iniset['specific_fermi']:
                self.fermi=iniset['specific_fermi_value']
            else:
                self.fermi=self.dos.fermi
        else:
            self.L_fermi=False
            self.fermi=0
        
    def plot_dos(self,plottype,pro_group=None,group_info=None,ax=None):
        if ax is None:
            ax=plt.subplot()
        if plottype==0:
            ax.plot(self.dos.Energy,self.dos.total.T-self.fermi)
        elif plottype==1:
            if (pro_group is None) | (group_info is None):
                raise ValueError(" There is no group info input")
            ax=self.plot_pdos(pro_group,group_info,ax=ax)
        if self.L_fermi:
            ax.plot([0,0],[-1000,1000],'--',c='gray')
        else:
            ax.plot([self.fermi,self.fermi],[-1000,1000],'--',c='gray')
        return ax

    def plot_pdos(self,pro_group,group_info,ax=None):
        if ax is None:
            ax=plt.subplot()
        for ispin in range(self.dos.I_spin):
            if ispin==0:
                label='Total'
            else:
                label=None
            ax.plot(self.dos.Energy-self.fermi,self.dos.total[ispin].T,c='gray')
            #plt.fill_between(self.dos.Energy-self.fermi,self.dos.total[ispin].T, color='gray', alpha=.25,label=label)
            for i,igroup in enumerate(pro_group):
                if ispin==0:
                    label=r''+group_info[i][0]
                else:
                    label=None
                ax.plot(self.dos.Energy-self.fermi,igroup[ispin].T,c=group_info[i][1])
                plt.fill_between(self.dos.Energy-self.fermi,igroup[ispin].T, color=group_info[i][1], alpha=.25,label=label)
        ax.legend(loc='right', bbox_to_anchor=(1, 0.6))
        return ax

class plot_chg(object):
    def __init__(self,chg):
        """
        this method only draw the image of one spin diagram
        """
        self.chg=chg
        pass
    def line_average(self,ispin,abc_axe,ax=plt.subplot()):
        """
        average alone a axis
        _axe = 
        'a','b','c' alone the lattice constant
        //'x','y','z' alone the x,y,z direction
        """
        _axe=['c','b','a'].index(abc_axe)
        #_axe=int(sys.argv[1])-1
        axe=[0,1,2]
        axe.remove(_axe)
        axe=tuple(axe)
        NG=self.chg.shape
        ax.plot(range(NG[_axe])/NG[_axe]*np.linalg.norm(self.chg.cell[(2,1,0)[_axe],:]),np.average(self.chg.chg,axis=axe),color='black')
        ax.set_xlim([np.min(range(NG[_axe])/NG[_axe]*np.linalg.norm(self.chg.cell[(2,1,0)[_axe],:])),np.max(range(NG[_axe])/NG[_axe]*np.linalg.norm(self.chg.cell[(2,1,0)[_axe],:]))])

