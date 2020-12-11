#!/usr/bin/python3
import sys

import matplotlib.pyplot as plt
import numpy as np


class get_chgcar(object):
    def __init__(self,LOCPOT='CHGCAR'):
        self.__file_name__=LOCPOT
        self.__file__=open(self.__file_name__)
        get_chgcar()

    def get_chgcar(self):
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        self.scale=float(tmp.split()[0])
        cell=np.zeros((3,3))
        for i in range(3):
            tmp=self.__file__.readline()
            cell[i,:]=tmp.split()[:3]
        self.element=self.__file__.readline().split()
        self.N_ions=np.array(self.__file__.readline().split(),dtype=np.int)
        tmp=self.__file__.readline()
        for i in range(np.sum(self.N_ions)):
            self.positions=np.array(self.__file__.readline().split(),dtype=np.float)
        tmp=self.__file__.readline()
        tmp=self.__file__.readline()
        self.NGX,self.NGY,self.NGZ=np.array(tmp.split()[:3],dtype='int')
        self.NG=(self.NGZ,self.NGY,self.NGX)
        LOC_DAT=[]
        for itmp in self.__file__:
            LOC_DAT.extend(itmp.split())
        try:
            self.chg=np.array(LOC_DAT,dtype='double').reshape((1,*self.NG))
        except:
            LOC_DATA=np.array(LOC_DAT[:np.prod(self.NG)],dtype='double').reshape(self.NG)
            LOC_DATB=np.array(LOC_DAT[-np.prod(self.NG):],dtype='double').reshape(self.NG)
            self.chg=np.vstack(((LOC_DATA+LOC_DATB)/2,(LOC_DATA-LOC_DATB)/2))

    def get_locpot(self):
        self.__init__('./POTCAR')

    def get_pchg(self):
        self.__init__('./PCHGCAR')

    def writen_spin_chg(self):
        self.__file__.seek(0)
        LOC=self.__file__
        LOC_a=open(self.__file_name__+"_alpha",'w')
        LOC_b=open(self.__file_name__+"_beta",'w')
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        cell=np.zeros((3,3))
        for i in range(3):
            tmp=LOC.readline()
            LOC_a.write(tmp)
            LOC_b.write(tmp)
            cell[i,:]=tmp.split()[:3]
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        NIONS=np.array(tmp.split(),dtype='int')
        NIONS_TOT=np.sum(NIONS)
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        for i in range(NIONS_TOT):
            tmp=LOC.readline()
            LOC_a.write(tmp)
            LOC_b.write(tmp)
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        tmp=LOC.readline()
        LOC_a.write(tmp)
        LOC_b.write(tmp)
        NGX,NGY,NGZ=np.array(tmp.split()[:3],dtype='int')
        NG=(NGZ,NGY,NGX)
        LOC_DAT=[]
        for itmp in LOC:
            LOC_DAT.extend(itmp.split())
        LOC_DAT=[np.array(LOC_DAT[:np.prod(NG)],dtype='double'),np.array(LOC_DAT[-np.prod(NG):],dtype='double')]
        LOC_DAT_a=(LOC_DAT[0]+LOC_DAT[1])/2
        LOC_DAT_b=(LOC_DAT[0]-LOC_DAT[1])/2
        count=0
        for temp in LOC_DAT_a:
            count+=1
            print("%12.9f" % temp,end=" ",file=LOC_a)
            if(count%5==0):
                print("", file=LOC_a)
        LOC_a.close()
        for temp in LOC_DAT_b:
            count+=1
            print("%12.9f" % temp,end=" ",file=LOC_b)
            if(count%5==0):
                print("",file=LOC_b,flush=True)
        LOC_b.close()

