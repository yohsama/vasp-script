#!/usr/bin/python3
import re

import numpy as np


class get_kpoints:
    def __init__(self,KPOINTS='./KPOINTS'):
        self.file=open(KPOINTS)
        tmp=self.file.readline()
        tmp=self.file.readline()
        self.N_kpt=int(tmp.split()[0])
        if self.N_kpt==0:
            print("It's a Auto-mesh KPOINTS file, trying to read IBZKPT")
            self.file.close()
            try:
                self.file=open('IBZKPT')
                tmp=self.file.readline()
                tmp=self.file.readline()
                self.N_kpt=int(tmp.split()[0])
            except:
                raise ValueError("IBZKPT is not found, please check the IBZKPT file or use a rec/line-mode format KPOINTS")
        tmp=self.file.readline()
        if 'l' in tmp[0].lower():
            print("Line-Mode KPOINTS file found")
            self.get_linemode_kpoints()
        elif 'r' in tmp[0].lower():
            print("rec KPOINTS file found")
            self.get_rec_kpoints()
        

    def get_rec_kpoints(self):
        self.file.seek(0)
        self.sp_kpt_label=[]
        self.file.readline()
        self.file.readline()
        self.file.readline()
        self.kpoint=np.zeros((self.N_kpt,3),dtype=float)
        self.kwht=np.zeros((self.N_kpt),dtype=float)
        tmp=self.file.readlines()
        inum=0
        for itmp in tmp:
            if inum<self.N_kpt:
                if re.match(' *-?[.0-9]{1,} *-?[.0-9]{1,} *-?[.0-9]{1,}',itmp):
                    self.kpoint[inum]=(np.array(itmp.split()[0:3],dtype=np.double))
                    self.kwht=(np.array(itmp.split()[3],dtype=np.double))
                    if "#" in itmp:
                        self.sp_kpt_label.append(r""+re.split(r"[# \n]+",itmp)[-2]+"")
                    else:
                        self.sp_kpt_label.append(r"")
                    inum=inum+1
        self.sp_kpt_label=np.array(self.sp_kpt_label)
        return self.kpoint,self.sp_kpt_label


    def get_linemode_kpoints(self):
        self.file.seek(0)
        self.sp_kpt_label=[]
        self.file.readline()
        self.file.readline()
        self.file.readline()
        tmp_kpoint=[]
        self.kwht=np.ones((self.N_kpt),dtype=float)
        tmp=self.file.readlines()
        tmp=np.array([ jtmp for jtmp in tmp if re.match(' *-?[.0-9]{1,} *-?[.0-9]{1,} *-?[.0-9]{1,}',jtmp) ])
        tmp=tmp.reshape((-1,2))
        for itmp in tmp:
            itmp1=itmp[0]
            itmp2=itmp[1]
            tmp1=np.array(itmp1.split()[0:3],dtype=np.double)
            tmp2=np.array(itmp2.split()[0:3],dtype=np.double)
            print(tmp1,tmp2)
            kx=np.interp(np.arange(0,self.N_kpt),[0,self.N_kpt],[tmp1[0],tmp2[0]])
            ky=np.interp(np.arange(0,self.N_kpt),[0,self.N_kpt],[tmp1[1],tmp2[1]])
            kz=np.interp(np.arange(0,self.N_kpt),[0,self.N_kpt],[tmp1[2],tmp2[2]])
            kv=np.vstack((kx,ky,kz)).T
            tmp_kpoint.append(kv)
            if "#" in itmp1:
                self.sp_kpt_label.append(r""+re.split(r"[# \n]+",itmp1)[-2]+"")
            else:
                self.sp_kpt_label.append(r"")
            self.sp_kpt_label.extend([r""]*(self.N_kpt-2))
            if "#" in itmp1:
                self.sp_kpt_label.append(r""+re.split(r"[# \n]+",itmp2)[-2]+"")
            else:
                self.sp_kpt_label.append(r"")
        self.sp_kpt_label=np.array(self.sp_kpt_label)
        self.kpoint=np.array(tmp_kpoint,dtype=np.double).reshape(-1,3)
        return self.kpoint,self.sp_kpt_label

    def make_k_path(self,kpoint,sp_kpt_label,bcell):
        kdiff=np.zeros(kpoint.shape)
        kdiff[1:]=np.diff(kpoint.dot(bcell),axis=0)
        kpath_section=np.linalg.norm(kdiff,axis=1)
        is_label=np.array(self.sp_kpt_label)!=""
        for i,tmp in enumerate((is_label[1:]) & (is_label[1:]==is_label[:-1])):
            if tmp:
                kpath_section[i+1]=0
                if self.sp_kpt_label[i+1]!=self.sp_kpt_label[i]:
                    self.sp_kpt_label[i]=self.sp_kpt_label[i]+r"|"+self.sp_kpt_label[i+1]
                self.sp_kpt_label[i+1]=""
        self.kpath=kpath_section.cumsum()
        return self.kpath,self.sp_kpt_label
    
