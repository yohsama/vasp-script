#!/usr/bin/env python
import numpy as np
from itertools import islice

class get_eigenvalue:
    def __init__(self,file="EIGENVAL"):
        self._file_=open(file)
        self.get_all_eigenvalue()
        pass

    def get_all_eigenvalue(self,occupied_threshold=0.001):
        self.N_ions=int(self._file_.readline().split()[0])
        tmp=self._file_.readline()
        tmp=self._file_.readline()
        tmp=self._file_.readline()
        tmp=self._file_.readline()
        tmp=self._file_.readline()
        self.N_Band=int(tmp.split()[2])
        self.N_kpt=int(tmp.split()[1])
        self.Kpoint=np.zeros((self.N_kpt,3),dtype=float)
        self.Kwht=np.zeros((self.N_kpt),dtype=float)
        EIG=[]
        for ikpt in range(self.N_kpt):
            self._file_.readline()
            self.Kpoint[ikpt,0],self.Kpoint[ikpt,1],self.Kpoint[ikpt,2],self.Kwht=np.array(self._file_.readline().split(),dtype=float)
            for itmp in islice(self._file_,0,self.N_Band):
                EIG.append(itmp.split())
        EIG=np.asarray(EIG,dtype=float).reshape(self.N_kpt,self.N_Band,-1)
        if EIG.shape[2] == 3:
            self.occ=EIG[:,:,2].reshape((1,self.N_kpt,-1))
            self.eig=EIG[:,:,1].reshape((1,self.N_kpt,-1))
            self.I_spin=1
        else:
            self.occ=EIG[:,:,3:].reshape((self.N_kpt,-1,2)).transpose((2,0,1))
            self.eig=EIG[:,:,1:3].reshape((self.N_kpt,-1,2)).transpose((2,0,1))
            self.I_spin=2

        self.fermi=np.max(self.eig[self.occ>occupied_threshold])
        return self.eig,self.occ,self.fermi
