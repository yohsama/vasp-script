#!/usr/bin/python3
import re
from itertools import islice

import numpy as np
import utilities


class get_procar:
    def __init__(self, PRO_FILES=['PROCAR']):
        self.__file__ = PRO_FILES
        for i, pfile in enumerate(PRO_FILES):
            if i == 0:
                self.eig, self.occ, self.project, self.kpoint, self.Kwht, self.N_ions = get_all_project_band(
                    pfile)
                self.N_spin, self.N_kpt, self.N_Band = self.eig.shape
                self.N_kpt = self.eig.shape
            else:
                eig, occ, kpoint, Kwht, N_ions = get_all_project_band(pfile)
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
                    self.project = np.append(self.occ, occ, axis=1)
                    self.kpoint = np.append(self.kpoint, kpoint, axis=0)
                    self.Kwht = np.append(self.Kwht, Kwht, axis=0)
                    self.N_kpt += eig.shape[1]

        self.fermi = utilities.get_fermi(self.eig, self.occ)

        return


def __get_L_orbit__(PROCAR='PROCAR'):
    with open(PROCAR) as __file__:
        tmp = __file__.readline()
        if "phase" in tmp:
            L_orbit = 12
            print("read as LORBIT == 12")
        elif "decomposed" in tmp:
            L_orbit = 11
            print("read as LORBIT == 11")
        elif "new" in tmp:
            L_orbit = 10
            print("read as LORBIT == 10")
        tmp = islice(__file__, 4, 5).__next__()
        if np.array(re.split(r'occ\.', tmp)[-1], dtype=np.float) > 1:
            I_spin = 1
        else:
            I_spin = 2
        if L_orbit > 10:
            return L_orbit, 9, I_spin
        elif L_orbit == 10:
            return L_orbit, 3, I_spin


def get_all_project_band(PROCAR='PROCAR'):
    L_orbit, __N_orbit__, I_spin = __get_L_orbit__(PROCAR)
    with open(PROCAR) as __file__:
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
                tmp = __file__.readline()
                kpoint[ikpt] = np.array([tmp[19:29], tmp[30:40], tmp[41:51]],
                                        dtype=np.float64)
                wht[ikpt] = np.float64(tmp[41:51])
                tmp = __file__.readline()  # 2
                for iband in range(N_band):
                    #print("ispin:%d,ikpt:%d,iband:%d" % (ispin,ikpt,iband))
                    tmp = __file__.readline()  #1*band
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
    return eig, occ, project, kpoint, wht, N_ions


def get_sp_project_band(PROCAR='PROCAR'):
    L_orbit, __N_orbit__, I_spin = __get_L_orbit__(PROCAR)
    with open(PROCAR) as __file__:
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
    return eig, occ, project, kpoint, wht, N_ions


def set_group(self, grouptag, symbollist, norm=True):
    return utilities.set_group(self, grouptag, symbollist, norm=True)
