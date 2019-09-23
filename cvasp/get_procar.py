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

def get_procar(PROCAR_files,ISPIN):
    i=0
    for PROCAR in PROCAR_files:
        f=open(PROCAR)
        tmp=f.readline()
        LORBIT=11
        if 'phase' in tmp:
            LORBIT=12
        tmp=f.readline()
        N_kpt=int(re.split(r'[^.0-9]+',tmp)[1])
        N_band=int(re.split(r'[^.0-9]+',tmp)[2])
        N_ions=int(re.split(r'[^.0-9]+',tmp)[3])
        tmp=f.readline()
        PRO_=np.zeros((ISPIN,N_kpt,N_band,N_ions))
        EIG_=np.zeros((ISPIN,N_kpt,N_band))
        OCC_=np.zeros((ISPIN,N_kpt,N_band))
        KPOINTS_=np.zeros((N_kpt,3))
        for ispin in range(ISPIN):
            for ikpt in range(N_kpt):
                tmp=f.readline()
                KPOINTS_[ikpt]=np.matrix([np.double(tmp[19:29]),np.double(tmp[30:40]),np.double(tmp[41:51])])
                tmp=f.readline()
                for iband in range(N_band):
                    tmp=f.readline()
                    EIG_[ispin,ikpt,iband]=re.split(r'[ \n]+',tmp)[4]
                    OCC_[ispin,ikpt,iband]=re.split(r'occ\.',tmp)[-1]
                    tmp=f.readline()
                    tmp=f.readline()
                    for iion in range(N_ions):
                        tmp=f.readline()
                        PRO_[ispin,ikpt,iband,iion]=tmp.split()[-1]
                    tmp=f.readline()
                    if LORBIT==12:
                        tmp=f.readline()
                        while not tmp.split()==[]:
                            tmp=f.readline()
                        #tmp=f.readline()
                    #tmp=f.readline()
                tmp=f.readline()
            f.readline()
        f.close()
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
    return EIG,OCC,PRO,KPOINTS
