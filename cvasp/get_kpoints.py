#!/usr/bin/python3
import numpy as np
import re


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
    print(K_num)
    sp_kpt=np.zeros((K_num,4))
    tmp=K.readlines();
    inum=0
    for itmp in tmp:
        if inum<K_num:
            if "." in itmp:
                sp_kpt[inum]=(np.array(itmp.split()[0:4],dtype="double"))
                if "#" in itmp:
                    sp_kpt_label.append(re.sub('\n','',re.split(r"[# ]",itmp)[-1]))
                else:
                    sp_kpt_label.append("")
                inum=inum+1
    K.close()
    return sp_kpt,sp_kpt_label


def make_k_path(KPOINTS_files,rep):
    kpt_label=[""]
    kpath=np.zeros((1,3))
    KPOINTS=np.zeros((1,3))
    for KPOINTS_file in KPOINTS_files:
        sp_kpt,sp_kpt_label=get_kpoints(KPOINTS_file)
        kpath=np.vstack((kpath,np.zeros((1,3)),sp_kpt[1:,0:3]-sp_kpt[:-1,0:3]))
        KPOINTS=np.vstack((KPOINTS,sp_kpt[:,0:3]))
        if (sp_kpt_label[0]!=kpt_label[-1]) & ((kpt_label[-1]!="")|(sp_kpt_label[0]!="")) & (len(kpt_label)>1):
            kpt_label[-1]=kpt_label[-1]+"|"+sp_kpt_label[0]
            sp_kpt_label[0]=kpt_label[-1]
            kpt_label.extend(sp_kpt_label)
        else:
            kpt_label.extend(sp_kpt_label)
    kpath=np.cumsum(np.linalg.norm(kpath.dot(rep),axis=1))
    kpath=kpath[1:]
    kpt_label=np.array(kpt_label)
    KPOINTS=KPOINTS[1:]
    kpt_label=kpt_label[1:]
    return kpath,kpt_label,KPOINTS
    
