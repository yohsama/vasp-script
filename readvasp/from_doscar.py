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

class get_doscar:
    def __init__(self,DOSCAR='./DOSCAR'):
        self.__file__=open(DOSCAR)
        self.get_dos()
    
    def get_dos(self):
        tmp=self.__file__.readline()
        self.N_ions=int(tmp.split()[0])
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        self.N_edos=int(tmp.split()[2])
        self.fermi=float(tmp.split()[3])
        tmp_total=[]
        for _ in range(self.N_edos):
            tmp=self.__file__.readline()
            tmp_total.append(np.array(tmp.split(),dtype=np.float))
        total=np.array(tmp_total,dtype=np.float).reshape(self.N_edos,-1)
        self.Energy=total[:,0]
        if total.shape[1]==3:
            self.I_spin=1
            self.total=total[:,1].reshape(1,-1)
            self.total_cumsum=total[:,2].reshape(1,-1)
        elif total.shape[1]==5:
            self.I_spin=2
            self.total=total[:,[1,2]].reshape(-1,2).transpose((1,0))
            self.total_cumsum=total[:,[3,4]].reshape(-1,2).transpose((1,0))
        else:
            print("Unrecognize Version in SPIN")
        self.project=[]
        for _ in range(self.N_ions):
            tmp=self.__file__.readline()    
            for _ in range(self.N_edos):
                tmp=self.__file__.readline()
                self.project.append(tmp.split()[1:])

        self.project=np.array(self.project,dtype=float)
        self.project=self.project.reshape((self.N_ions,self.N_edos,-1,self.I_spin))
        self.project=self.project.transpose((3,1,0,2))
        if self.project.shape[3] == 3:
            self.L_orbit=10
        elif self.project.shape[3] == 9:
            self.L_orbit=11
        else:
            print("Unrecognize Version in LORBIT; code [%d]" % self.project.shape[2])
        if self.I_spin==2:
            self.total_cumsum[1]=-self.total_cumsum[1]
            self.total[1]=-self.total[1]
            self.project[1]=-self.project[1]

    def set_group(self,grouptag,symbollist,norm=True):
        return ultils.set_group(self,grouptag,symbollist,norm=True)