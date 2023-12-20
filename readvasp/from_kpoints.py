#!/usr/bin/python3
import re

import ase.io
import numpy as np


class get_kpoints(object):
    def __init__(self, KPOINTS_FILE=['KPOINTS']):
        self.ini_kpt()
        self.kfiles = KPOINTS_FILE
        for kfile in KPOINTS_FILE:
            print('reading %s' % kfile)
            self.gen_kpt(kfile)
        # self.k_sp_label=np.array(self.k_sp_label)
        self.make_k_path()

    def ini_kpt(self):
        self.k_sp_label = []
        #self.kpoint = np.zeros((0, 3), dtype=float)
        self.N_kpt = 0
        self.kpoint=np.empty((0,3))
        self.kwht=np.empty(0)

    def gen_kpt(self, KPOINTS='./KPOINTS'):
        self.file = open(KPOINTS)
        tmp = self.file.readline()
        tmp = self.file.readline()
        if int(tmp.split()[0]) == 0:
            print("It's a Auto-mesh KPOINTS file, trying to read IBZKPT")
            self.file.close()
            try:
                self.file = open('IBZKPT')
                tmp = self.file.readline()
                tmp = self.file.readline()
            except:
                raise ValueError(
                    "IBZKPT is not found, please check the IBZKPT file or use a rec/line-mode format KPOINTS"
                )
        self.N_kpt += int(tmp.split()[0])
        tmp = self.file.readline()
        if 'l' in tmp[0].lower():
            print("Line-Mode KPOINTS file found")
            self.get_linemode_kpoints()
        elif 'r' in tmp[0].lower():
            print("rec KPOINTS file found")
            self.get_rec_kpoints()

    def get_rec_kpoints(self):
        self.file.seek(0)
        self.file.readline()
        tmp=self.file.readline()
        N_kpt=int(tmp.split()[0])
        self.file.readline()
        kwht = np.zeros((N_kpt), dtype=float)
        kpoint = np.zeros((N_kpt, 3), dtype=float)
        tmp = self.file.readlines()
        inum = 0
        for itmp in tmp:
            if ((inum < N_kpt) & bool(re.match(' *-?[.0-9]{1,} *-?[.0-9]{1,} *-?[.0-9]{1,}', itmp))):
                kpoint[inum] = np.array(itmp.split()[:3], dtype=np.float64)
                kwht[inum] = (np.array(itmp.split()[3], dtype=np.float64))
                if "#" in itmp:
                    self.k_sp_label.append(
                        r"" + re.split(r"[# \n]+", itmp)[-2] + "")
                else:
                    self.k_sp_label.append(r"")
                inum = inum + 1
        self.kpoint=np.concatenate((self.kpoint,kpoint))
        self.kwht=np.concatenate((self.kwht,kwht))
    
    def get_linemode_kpoints(self):
        self.file.seek(0)
        self.file.readline()
        self.file.readline()
        self.file.readline()
        kpoint = []
        kwht = np.ones((self.N_kpt), dtype=float)
        tmp = self.file.readlines()
        tmp = np.array([
            jtmp for jtmp in tmp
            if re.match(' *-?[.0-9]{1,} *-?[.0-9]{1,} *-?[.0-9]{1,}', jtmp)
        ])
        tmp = tmp.reshape((-1, 2))
        for itmp in tmp:
            itmp1 = itmp[0]
            itmp2 = itmp[1]
            tmp1 = np.array(itmp1.split()[0:3], dtype=np.float64)
            tmp2 = np.array(itmp2.split()[0:3], dtype=np.float64)
            #print(tmp1, tmp2)
            kx = np.interp(np.arange(0, self.N_kpt), [0, self.N_kpt],
                           [tmp1[0], tmp2[0]])
            ky = np.interp(np.arange(0, self.N_kpt), [0, self.N_kpt],
                           [tmp1[1], tmp2[1]])
            kz = np.interp(np.arange(0, self.N_kpt), [0, self.N_kpt],
                           [tmp1[2], tmp2[2]])
            kv = np.vstack((kx, ky, kz)).T
            kpoint.append(kv)
            if "#" in itmp1:
                self.k_sp_label.append(r"" +
                                       re.split(r"[# \n]+", itmp1)[-2] + "")
            else:
                self.k_sp_label.append(r"")
            self.k_sp_label.extend([r""] * (self.N_kpt - 2))
            if "#" in itmp2:
                self.k_sp_label.append(r"" +
                                       re.split(r"[# \n]+", itmp2)[-2] + "")
            else:
                self.k_sp_label.append(r"")
        print(self.kpoint.shape,np.array(kpoint).shape)
        self.kpoint=np.concatenate((self.kpoint,np.array(kpoint).reshape(-1,3)))
        self.kwht=np.concatenate((self.kwht,kwht))
        
    def make_k_path(self):
        root = re.sub('KPOINTS(^/)*$', '', self.kfiles[0])
        for str_file in ['CONTCAR', 'POSCAR', 'vasprun.xml', 'OUTCAR']:
            try:
                pos = ase.io.read(root+"/"+str_file, "-1")
                bcell = pos.cell.reciprocal()
            except:
                pass
            else:
                print("structure read from %s" % str_file)
                break
        else:
            print("No structure file found")
            bcell = np.identity(3)
        self.kpath, self.k_sp_label = make_k_path(self.kpoint, self.k_sp_label,
                                                  bcell)


def make_k_path(kpoint, k_sp_label, bcell):
    kdiff = np.zeros(kpoint.shape)
    kdiff[1:] = np.diff(kpoint.dot(bcell), axis=0,)
    kpath_section = np.linalg.norm(kdiff, axis=1)
    for i in range(len(k_sp_label) - 1):
        if bool(k_sp_label[i]) & (not "!" in k_sp_label[i]):
            # print(k_sp_label[i])
            if k_sp_label[i] == k_sp_label[i + 1]:
                k_sp_label[i] = ""
                kpath_section[i + 1] = 0
            elif bool(k_sp_label[i+1]) & (not "!" in k_sp_label[i+1]):
                kpath_section[i + 1] = 0
                k_sp_label[i] = k_sp_label[i] + r"|" + k_sp_label[i + 1]
                k_sp_label[i + 1]=""
    kpath = kpath_section.cumsum()
    return kpath, k_sp_label
