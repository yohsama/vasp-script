#!/usr/bin/env python
import numpy as np
from itertools import islice


class get_eigenvalue:
    def __init__(self, EIG_FILES=["EIGENVAL"]):
        self.__file__ = EIG_FILES
        for i, efile in enumerate(EIG_FILES):
            print('reading %s' % efile)
            if i == 0:
                self.eig, self.occ, self.kpoint, self.Kwht, self.N_ions = __get_eigenvalue__(
                    efile)
                self.N_spin, self.N_kpt, self.N_Band = self.eig.shape
            else:
                eig, occ, kpoint, Kwht, N_ions = __get_eigenvalue__(efile)
                if (self.N_ions != N_ions):
                    print("Different ions found, last file is %s" % efile)
                elif (self.N_Band != eig.shape[2]):
                    print("Different band number found, last file is %s" %
                          efile)
                elif (self.N_spin != eig.shape[0]):
                    print("Different spin found, last file is %s" % efile)
                else:
                    self.eig = np.append(self.eig, eig, axis=1)
                    self.occ = np.append(self.occ, occ, axis=1)
                    self.kpoint = np.append(self.kpoint, kpoint, axis=0)
                    #self.Kwht = np.append(self.Kwht, Kwht, axis=0)
                    self.N_kpt += eig.shape[1]
            print(self.Kwht)
        self.fermi = get_fermi(self.eig, self.occ)


def __get_eigenvalue__(file='EIGENVAL'):
    with open(file) as __file__:
        N_ions = int(__file__.readline().split()[0])
        tmp = __file__.readline()
        tmp = __file__.readline()
        tmp = __file__.readline()
        tmp = __file__.readline()
        tmp = __file__.readline()
        N_Band = int(tmp.split()[2])
        N_kpt = int(tmp.split()[1])
        kpoint = np.zeros((N_kpt, 3), dtype=float)
        Kwht = np.zeros((N_kpt), dtype=float)
        EIG = []
        for ikpt in range(N_kpt):
            __file__.readline()
            kpoint[ikpt, 0], kpoint[ikpt, 1], kpoint[ikpt, 2], Kwht = np.array(
                __file__.readline().split(), dtype=float)
            for itmp in islice(__file__, 0, N_Band):
                EIG.append(itmp.split())
        EIG = np.asarray(EIG, dtype=float).reshape(N_kpt, N_Band, -1)
        if EIG.shape[2] == 3:
            occ = EIG[:, :, 2].reshape((1, N_kpt, -1))
            eig = EIG[:, :, 1].reshape((1, N_kpt, -1))
            print("find ispin = 1")
        else:
            occ = EIG[:, :, 3:].reshape((N_kpt, -1, 2)).transpose((2, 0, 1))
            eig = EIG[:, :, 1:3].reshape((N_kpt, -1, 2)).transpose((2, 0, 1))
            print("find ispin = 2")
    return eig, occ, kpoint, Kwht, N_ions


def get_fermi(eig, occ, occupied_threshold=0.001):
    return np.max(eig[occ > occupied_threshold])
