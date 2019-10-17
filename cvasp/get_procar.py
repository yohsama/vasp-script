#!/usr/bin/python3
import ase.io
import numpy as np
import re
import matplotlib.pyplot as plt
from glob  import glob
from scipy.interpolate import griddata
import os
import sys
import pandas as pd

def get_procar(PROCAR_files,ISPIN):
    i=0
 #   procar=pd.DataFrame(columns=('ispin','kpt','band','ion','occ','eig','proj_tot','proj_s','proj_py','proj_pz','proj_px','proj_dxy','proj_dyz','proj_dz2','proj_dxy','proj_dx2-y2');
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
        PRO_=np.zeros((ISPIN,N_kpt,N_band,N_ions,10))
        EIG_=np.zeros((ISPIN,N_kpt,N_band))
        OCC_=np.zeros((ISPIN,N_kpt,N_band))
        TTOT_=np.zeros((ISPIN,N_kpt,N_band))
        KPOINTS_=np.zeros((N_kpt,3))
        for ispin in range(ISPIN):
            for ikpt in range(N_kpt):
                tmp=f.readline()
                KPOINTS_[ikpt]=np.matrix([np.double(tmp[19:29]),np.double(tmp[30:40]),np.double(tmp[41:51])])
                tmp=f.readline()
                for iband in range(N_band):
                    tmp=f.readline()
                    #print(re.split(r'[ \n]+',tmp))
                    EIG_[ispin,ikpt,iband]=re.split(r'[ \n]+',tmp)[4]
                    OCC_[ispin,ikpt,iband]=re.split(r'occ\.',tmp)[-1]
                    tmp=f.readline()
                    tmp=f.readline()
                    for iion in range(N_ions):
                        tmp=f.readline()
                        PRO_[ispin,ikpt,iband,iion,0:len(tmp.split())-1]=tmp.split()[1:]
                        PRO_[ispin,ikpt,iband,iion,9]=tmp.split()[-1]
                    tmp=f.readline()
                    TTOT_[ispin,ikpt,iband]=tmp.split()[-1]
                    if LORBIT==12:
                        tmp=f.readline()
                        while not tmp.split()==[]:
                            tmp=f.readline()
                    else:
                        f.readline()
                tmp=f.readline()
            f.readline()
        f.close()
        if i==0:
            EIG=EIG_
            PRO=PRO_
            OCC=OCC_
            KPOINTS=KPOINTS_
            TTOT=TTOT_
        else:
            EIG=np.concatenate((EIG,EIG_),axis=1)
            PRO=np.concatenate((PRO,PRO_),axis=1)
            OCC=np.concatenate((OCC,OCC_),axis=1)
            TTOT=np.concatenate((TTOT,TTOT_),axis=1)
            KPOINTS=np.concatenate((KPOINTS,KPOINTS_),axis=0)
        i=i+1
        spin    = np.arange(PRO.shape[0]).repeat(PRO.shape[1]).repeat(PRO.shape[2]).repeat(PRO.shape[3]).flatten()
        kpt     = np.tile(np.arange(PRO.shape[1]).repeat(PRO.shape[2]).repeat(PRO.shape[3]),PRO.shape[0]).flatten()
        band    = np.tile(np.tile(np.arange(PRO.shape[2]).repeat(PRO.shape[3]),PRO.shape[1]),PRO.shape[0]).flatten()+1
        ion     = np.tile(np.tile(np.tile(np.arange(PRO.shape[3]),PRO.shape[2]),PRO.shape[1]),PRO.shape[0]).flatten()+1
        ttot    = TTOT.repeat(PRO.shape[3])
        KPOINTS1= np.tile(KPOINTS[:,0].repeat(N_ions).repeat(N_band),ISPIN).flatten()
        KPOINTS2= np.tile(KPOINTS[:,1].repeat(N_ions).repeat(N_band),ISPIN).flatten()
        KPOINTS3= np.tile(KPOINTS[:,2].repeat(N_ions).repeat(N_band),ISPIN).flatten()
        occ     = OCC.repeat(PRO.shape[3])
        eig     = EIG.repeat(PRO.shape[3])
        procar=pd.DataFrame(PRO.reshape((int(PRO.size/PRO.shape[-1]),PRO.shape[-1])),columns=('s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','tot')  )
        procar['ttot']=ttot
        procar['spin']=spin
        procar['kpt']=kpt
        procar['band']=band
        procar['ion']=ion
        procar['occ']=occ
        procar['eig']=eig
        procar['kpoints1']=KPOINTS1
        procar['kpoints2']=KPOINTS2
        procar['kpoints3']=KPOINTS3
    return procar
