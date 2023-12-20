#!/usr/bin/python3
import sys

import matplotlib.pyplot as plt
import numpy as np


class get_chgcar(object):
    def __init__(self, LOCPOT='CHGCAR'):
        self.__file_name__ = LOCPOT
        self.__file__ = open(self.__file_name__)
        self.get_chgcar()

    def get_chgcar(self):
        self.system = self.__file__.readline()[:-2]
        self._factor_ = float(self.__file__.readline().split()[0])
        self.cell = np.zeros((3, 3))
        for i in range(3):
            self.cell[i, :] = self.__file__.readline().split()[:3]
        self.element = self.__file__.readline().split()
        self.N_ions = np.array(self.__file__.readline().split(), dtype=np.int)
        self.type = self.__file__.readline()[:-2]
        self.positions = np.array([
            np.array(self.__file__.readline().split()[:3], dtype=np.float64)
            for i in range(np.sum(self.N_ions))
        ])

        self.__file__.readline()
        self.NG = np.array(self.__file__.readline().split()[:3], dtype='int')
        LOC_DAT = []
        for itmp in self.__file__:
            LOC_DAT.extend(itmp.split())
        try:
            self.chg = np.array(LOC_DAT[:np.prod(self.NG)],
                                dtype=np.float64).reshape((1, *self.NG),
                                                         order='F')
        except:
            LOC_DATA = np.array(LOC_DAT[:np.prod(self.NG)],
                                dtype=np.float64).reshape(self.NG, order='F')
            LOC_DATB = np.array(
                LOC_DAT[len(LOC_DAT) // 2 + 3:len(LOC_DAT) // 2 +
                        np.prod(self.NG) + 3],
                dtype=np.float64).reshape(self.NG, order='F')
            self.chg = np.stack(
                ((LOC_DATA + LOC_DATB) / 2, (LOC_DATA - LOC_DATB) / 2))

    def write_spin_chg(self):
        write_chgcar(self, self.chg[0], self.__file__.name + "_alpha")
        write_chgcar(self, self.chg[1], self.__file__.name + "_beta")


def write_chgcar(chg, data, file):
    with open(file, 'w') as f:
        print(chg.system, file=f)
        print(' {: 18.15f}'.format(chg._factor_), file=f)
        np.savetxt(f, chg.cell, fmt='  % 11.6f% 11.6f% 11.6f')
        print(*chg.element, file=f)
        np.savetxt(f, chg.N_ions[None, :], fmt='%d')
        print(chg.type, file=f)
        np.savetxt(f, chg.positions, fmt='% 10.6f')
        print("", file=f)
        np.savetxt(f, np.array(chg.NG)[None, :], fmt='%d')
        cols = data.size // 5
        np.savetxt(f,
                   data.reshape(-1, order='F')[:5 * cols].reshape(-1, 5),
                   fmt='% 18.11e')
        np.savetxt(f,
                   data.reshape(-1, order='F')[None, 5 * cols:],
                   fmt='% 18.11e')


def diff1(chg_total, chg_A, chg_B):
    data = chg_total.chg - chg_A.chg - chg_B.chg
    write_chgcar(chg_total, data, chg_total.__file__.name + "_diff")


def diff2(chg_total, chg_A, chg_B, chg_C):
    data = chg_total.chg - chg_A.chg - chg_B.chg - chg_C.chg
    write_chgcar(chg_total, data, chg_total.__file__.name + "_diff")
