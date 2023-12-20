#!/usr/bin/python3
import re
from itertools import islice

import numpy as np
import utilities
import gzip

class get_procar(object):
    def __init__(self, PRO_FILES=['PROCAR']):
        self.__file__ = PRO_FILES
        for i, pfile in enumerate(PRO_FILES):
            print('reading %s' % pfile)
            if i == 0:
                self.L_orbit, self.eig, self.occ, self.project, self.kpoint, self.Kwht, self.N_ions = self.get_all_project_band(
                    pfile)
                self.N_spin, self.N_kpt, self.N_Band = self.eig.shape
            else:
                L_orbit, eig, occ, project, kpoint, Kwht, N_ions = self.get_all_project_band(pfile)
                if (self.N_ions != N_ions):
                    print("Different ions found, last file is %s" % pfile)
                elif (self.N_Band != eig.shape[2]):
                    print("Different band number found, last file is %s" %
                          pfile)
                elif (self.N_spin != eig.shape[0]):
                    print("Different spin found, last file is %s" % pfile)
                else:
                    self.eig = np.append(self.eig, eig, axis=1)
                    self.occ = np.append(self.occ, occ, axis=1)
                    self.project = np.append(self.project, project, axis=1)
                    self.kpoint = np.append(self.kpoint, kpoint, axis=0)
                    self.Kwht = np.append(self.Kwht, Kwht, axis=0)
                    self.N_kpt += eig.shape[1]

        self.fermi = utilities.get_fermi(self.eig, self.occ)

    def set_group(self, grouptag, symbollist):
        return utilities.set_group(self.L_orbit, self.project, grouptag, symbollist)

    def __get_L_orbit__(self,PROCAR='PROCAR'):
        with gzip.open(PROCAR,mode='rt') if PROCAR.endswith('.gz') else open(PROCAR) as __file__:
            tmp = __file__.readline()
            print(tmp)
            if "phase" in tmp:
                self.L_orbit = 12
                print("read as LORBIT == 12")
            elif "decomposed" in tmp:
                self.L_orbit = 11
                print("read as LORBIT == 11")
            elif "new" in tmp:
                self.L_orbit = 10
                print("read as LORBIT == 10")
            tmp = islice(__file__, 4, 5).__next__()
            if np.array(re.split(r'occ\.', tmp)[-1], dtype=np.float) > 1:
                self.I_spin = 1
            else:
                self.I_spin = 2
            tmp = islice(__file__,1,2).__next__()
            
            self.__N_orbit__ = len(tmp.split())-2
            # if self.L_orbit > 10:
            #     self.__N_orbit__=9
            # elif self.L_orbit == 10:
            #     self.__N_orbit__=3
            # print(self.__N_orbit__)

    def get_all_project_band(self,PROCAR='PROCAR'):
        self.__get_L_orbit__(PROCAR)
        L_orbit, __N_orbit__, I_spin = self.L_orbit, self.__N_orbit__, self.I_spin 
        with gzip.open(PROCAR,mode='rt') if PROCAR.endswith('.gz') else open(PROCAR) as __file__:

                tmp = __file__.readline()
                tmp = __file__.readline()
                N_kpt = int(re.split(r'[^.0-9]+', tmp)[1])
                N_band = int(re.split(r'[^.0-9]+', tmp)[2])
                N_ions = int(re.split(r'[^.0-9]+', tmp)[3])
                project = np.zeros((I_spin, N_kpt, N_band, N_ions, __N_orbit__))
                total_of_project_atom = np.zeros((I_spin, N_kpt, N_band, N_ions))
                total_of_project_orbit = np.zeros((I_spin, N_kpt, N_band, __N_orbit__))
                eig = np.zeros((I_spin, N_kpt, N_band))
                occ = np.zeros((I_spin, N_kpt, N_band))
                total = np.zeros((I_spin, N_kpt, N_band))
                kpoint = np.zeros((N_kpt, 3))
                wht = np.zeros((N_kpt))
                tmp = __file__.readline()
                for ispin in range(I_spin):
                    for ikpt in range(N_kpt):
                        while not "k-point " in tmp:
                            tmp = __file__.readline() 
                        kpoint[ikpt] = np.array([tmp[19:29], tmp[30:40], tmp[41:51]],
                                                dtype=np.float64)
                        wht[ikpt] = np.float64(tmp[41:51])
                        # tmp = __file__.readline()  # 2
                        for iband in range(N_band):
                            #print("ispin:%d,ikpt:%d,iband:%d" % (ispin,ikpt,iband))
                            while not "band" in tmp:
                                tmp = __file__.readline() 
                            eig[ispin, ikpt, iband] = re.split(r'[ \n]+', tmp)[4]
                            occ[ispin, ikpt, iband] = re.split(r'occ\.', tmp)[-1]
                            tmp = __file__.readline()  #2*band
                            tmp = __file__.readline()  #3*band
                            for iion in range(N_ions):
                                tmp = __file__.readline()  # (3+N_ions)*band
                                project[ispin, ikpt, iband, iion] = tmp.split()[1:-1]
                                total_of_project_atom[ispin, ikpt, iband,
                                                    iion] = tmp.split()[-1]
                            if N_ions > 1:
                                tmp = __file__.readline()  # (4+N_ions)*band
                                total_of_project_orbit[ispin, ikpt,
                                                    iband] = tmp.split()[1:-1]
                                total[ispin, ikpt, iband] = tmp.split()[-1]
                            if L_orbit == 12:
                                tmp = __file__.readline()  #(4/5+N_ions)*band
                                lmdata = []
                                while not tmp.split() == []:
                                    tmp = __file__.readline()
                                    lmdata.append(tmp)  #(4/5+N_ions)*band
                            else:
                                __file__.readline()
                        tmp = __file__.readline()
                    __file__.readline()
                print("PROCAR read complete")
        return L_orbit, eig, occ, project, kpoint, wht, N_ions


    def get_sp_project_band(self, PROCAR='PROCAR'):
        self.__get_L_orbit__(PROCAR)
        L_orbit, __N_orbit__, I_spin = self.L_orbit, self.__N_orbit__, self.I_spin 
        with gzip.open(PROCAR,mode='rt') if PROCAR.endswith('.gz') else open(PROCAR) as __file__:
     
                tmp = __file__.readline()
                tmp = __file__.readline()
                N_kpt = int(re.split(r'[^.0-9]+', tmp)[1])
                N_band = int(re.split(r'[^.0-9]+', tmp)[2])
                N_ions = int(re.split(r'[^.0-9]+', tmp)[3])
                project = np.zeros((I_spin, N_kpt, N_band, N_ions, 10))
                eig = np.zeros((I_spin, N_kpt, N_band))
                occ = np.zeros((I_spin, N_kpt, N_band))
                total = np.zeros((I_spin, N_kpt, N_band))
                kpoint = np.zeros((N_kpt, 3))
                wht = np.zeros((N_kpt))
                tmp = __file__.readline()
                for ispin in range(I_spin):
                    for ikpt in range(N_kpt):
                        tmp = __file__.readline()
                        kpoint[ikpt] = np.array([tmp[19:29], tmp[30:40], tmp[41:51]],
                                                dtype=np.float64)
                        wht[ikpt] = np.float64(tmp[41:51])
                        tmp = __file__.readline()
                        for iband in range(N_band):
                            tmp = __file__.readline()
                            eig[ispin, ikpt, iband] = re.split(r'[ \n]+', tmp)[4]
                            occ[ispin, ikpt, iband] = re.split(r'occ\.', tmp)[-1]
                            tmp = __file__.readline()
                            tmp = __file__.readline()
                            for iion in range(N_ions):
                                tmp = __file__.readline()
                                project[ispin, ikpt, iband, iion,
                                        0:len(tmp.split()) - 1] = tmp.split()[1:]
                                project[ispin, ikpt, iband, iion, -1] = tmp.split()[-1]
                            if N_ions > 1:
                                tmp = __file__.readline()
                                total[ispin, ikpt, iband] = tmp.split()[-1]
                            if L_orbit == 12:
                                tmp = __file__.readline()
                                while not tmp.split() == []:
                                    tmp = __file__.readline()
                            else:
                                __file__.readline()
                        tmp = __file__.readline()
                    __file__.readline()
                print("PROCAR read complete")
        return L_orbit, eig, occ, project, kpoint, wht, N_ions



