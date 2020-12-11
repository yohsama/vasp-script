#!/usr/bin/python3
import os
import re
import sys
from glob import glob
from itertools import islice

import ase.io
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import  ultils


class get_procar:
    def __init__(self,PROCAR='./PROCAR',occupied_threshold=0.001):
        self.file=open(PROCAR)
        self.file.seek(0)
        tmp=islice(self.file,5,6).__next__()
        if np.array(re.split(r'occ\.',tmp)[-1],dtype=np.float ) > 1:
            self.I_spin = 1
        else:
            self.I_spin = 2
        self.get_L_orbit()
        self.get_all_project_band()
        self.fermi=np.max(self.eig[self.occ>occupied_threshold])
        print("PROCAR read complete")
        pass

    def get_L_orbit(self):
        self.file.seek(0)
        tmp=self.file.readline()
        if "phase" in tmp:
            self.L_orbit=12
            print("read as LORBIT == 12")
        elif "decomposed" in tmp:
            self.L_orbit=11
            print("read as LORBIT == 11")
        elif "new" in tmp:
            self.L_orbit=10
            print("read as LORBIT == 10")

    def get_all_project_band(self):
        self.file.seek(0)
        tmp=self.file.readline()
        tmp=self.file.readline()
        self.N_kpt=int(re.split(r'[^.0-9]+',tmp)[1])
        self.N_band=int(re.split(r'[^.0-9]+',tmp)[2])
        self.N_ions=int(re.split(r'[^.0-9]+',tmp)[3])
        if self.L_orbit>10:
            self.__N_orbit__=9
        elif self.L_orbit==10:
            self.__N_orbit__=3
        self.project=np.zeros((self.I_spin,self.N_kpt,self.N_band,self.N_ions,self.__N_orbit__))
        self.total_of_project_atom=np.zeros((self.I_spin,self.N_kpt,self.N_band,self.N_ions))
        self.total_of_project_orbit=np.zeros((self.I_spin,self.N_kpt,self.N_band,self.__N_orbit__))
        self.eig=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.occ=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.total=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.kpoint=np.zeros((self.N_kpt,3))
        self.wht=np.zeros((self.N_kpt))
        tmp=self.file.readline() 
        for ispin in range(self.I_spin):
            for ikpt in range(self.N_kpt):
                tmp=self.file.readline() 
                self.kpoint[ikpt]=np.array([tmp[19:29],tmp[30:40],tmp[41:51]],dtype=np.double)
                self.wht[ikpt]=np.double(tmp[41:51])
                tmp=self.file.readline() # 2
                for iband in range(self.N_band):
                    #print("ispin:%d,ikpt:%d,iband:%d" % (ispin,ikpt,iband))
                    tmp=self.file.readline() #1*band 
                    self.eig[ispin,ikpt,iband]=re.split(r'[ \n]+',tmp)[4]
                    self.occ[ispin,ikpt,iband]=re.split(r'occ\.',tmp)[-1]
                    tmp=self.file.readline() #2*band
                    tmp=self.file.readline() #3*band
                    for iion in range(self.N_ions):
                        tmp=self.file.readline() # (3+self.N_ions)*band
                        self.project[ispin,ikpt,iband,iion]=tmp.split()[1:-1]
                        self.total_of_project_atom[ispin,ikpt,iband,iion]=tmp.split()[-1]
                    if self.N_ions>1:
                        tmp=self.file.readline() # (4+self.N_ions)*band
                        self.total_of_project_orbit[ispin,ikpt,iband]=tmp.split()[1:-1]
                        self.total[ispin,ikpt,iband]=tmp.split()[-1]
                    if self.L_orbit==12:
                        tmp=self.file.readline() #(4/5+self.N_ions)*band
                        self.lmdata=[]
                        while not tmp.split()==[]:
                            tmp=self.file.readline()
                            self.lmdata.append(tmp) #(4/5+self.N_ions)*band
                    else:
                        self.file.readline()
                tmp=self.file.readline()
            self.file.readline()
        return
    
    def get_sp_project_band(self):
        self.file.seek(0)
        tmp=self.file.readline()
        tmp=self.file.readline()
        self.N_kpt=int(re.split(r'[^.0-9]+',tmp)[1])
        self.N_band=int(re.split(r'[^.0-9]+',tmp)[2])
        self.N_ions=int(re.split(r'[^.0-9]+',tmp)[3])
        self.project=np.zeros((self.I_spin,self.N_kpt,self.N_band,self.N_ions,10))
        self.eig=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.occ=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.total=np.zeros((self.I_spin,self.N_kpt,self.N_band))
        self.kpoint=np.zeros((self.N_kpt,3))
        self.wht=np.zeros((self.N_kpt))
        tmp=self.file.readline() 
        for ispin in range(self.I_spin):
            for ikpt in range(self.N_kpt):
                tmp=self.file.readline()
                self.kpoint[ikpt]=np.array([tmp[19:29],tmp[30:40],tmp[41:51]],dtype=np.double)
                self.wht[ikpt]=np.double(tmp[41:51])
                tmp=self.file.readline()
                for iband in range(self.N_band):
                    tmp=self.file.readline()
                    self.eig[ispin,ikpt,iband]=re.split(r'[ \n]+',tmp)[4]
                    self.occ[ispin,ikpt,iband]=re.split(r'occ\.',tmp)[-1]
                    tmp=self.file.readline()
                    tmp=self.file.readline()
                    for iion in range(self.N_ions):
                        tmp=self.file.readline()
                        self.project[ispin,ikpt,iband,iion,0:len(tmp.split())-1]=tmp.split()[1:]
                        self.project[ispin,ikpt,iband,iion,-1]=tmp.split()[-1]
                    if self.N_ions>1:
                        tmp=self.file.readline()
                        self.total[ispin,ikpt,iband]=tmp.split()[-1]
                    if self.L_orbit==12:
                        tmp=self.file.readline()
                        while not tmp.split()==[]:
                            tmp=self.file.readline()
                    else:
                        self.file.readline()
                tmp=self.file.readline()
            self.file.readline()
        return
    
    def set_group(self,grouptag,symbollist,norm=True):
        return ultils.set_group(self,grouptag,symbollist,norm=True)

