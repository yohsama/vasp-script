#!/usr/bin/python3
import numpy as np
import re
import pandas as pd

def get_kpoints(KPOINTS):
    K=open(KPOINTS)
    sp_kpt_label=[]
    K.readline()
    tmp=K.readline()
    K_num=int(tmp.split()[0])
    if K_num==0:
        K.close()
        K=open('IBZKPT')
        K.readline()
        tmp=K.readline()
        K_num=int(tmp.split()[0])
    sp_kpt=np.zeros((K_num,4))
    tmp=K.readlines();
    inum=0
    for itmp in tmp:
        if inum<K_num:
            if "." in itmp:
                sp_kpt[inum]=(np.array(itmp.split()[0:4],dtype="double"))
                if "#" in itmp:
                    sp_kpt_label.append(r""+re.split(r"[# \n]+",itmp)[-2]+"")
                else:
                    sp_kpt_label.append(r"")
                inum=inum+1
    K.close()
    return sp_kpt,sp_kpt_label


def make_k_path(KPOINTS_files,rep):
    kpt_label=[""]
    kpath=np.zeros((1,3))
    KPOINTS=np.zeros((1,3))
    lastlabel=""
    last=0
    for KPOINTS_file in KPOINTS_files:
        sp_kpt,sp_kpt_label=get_kpoints(KPOINTS_file)
        kpath=np.vstack((kpath,np.zeros((1,3)),sp_kpt[1:,0:3]-sp_kpt[:-1,0:3]))
        KPOINTS=np.vstack((KPOINTS,sp_kpt[:,0:3]))
        kpt_label.extend(sp_kpt_label)
    for i in range(len(kpt_label)):
        if "!" in kpt_label[i]:
            try:
                kpath[i+1]=0
            except:
                pass
            kpath[i]=0
        else :
            if  (kpt_label[i]!=lastlabel) & ((lastlabel!="")&(kpt_label[i]!="")) & (len(kpt_label)>1) :
                kpt_label[last]=""
                kpt_label[i]=lastlabel+r'$\vert$'+kpt_label[i]
            lastlabel=kpt_label[i]
            last=i
    kpath=np.cumsum(np.linalg.norm(kpath.dot(rep),axis=1))
    kpath=kpath[1:]
    kpt_label=np.array(kpt_label)
    KPOINTS=KPOINTS[1:]
    kpt_label=kpt_label[1:]
    kpoints=pd.DataFrame(KPOINTS,columns=('kpoint1','kpoint2','kpoint3'))
    kpoints['kpath']=kpath
    kpoints['kpt_label']=kpt_label
    kpoints['kpt']=np.arange(kpath.size)
    return kpoints
    
