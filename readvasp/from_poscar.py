
import ase.io
import numpy as np


class get_poscar(object):
    def __init__(self,POSCAR='./POSCAR'):
        self.__file_name__=POSCAR
        self.get_by_ase()
        pass

    def get_symbollist(self):
        symbol=np.array(self.atoms.get_chemical_symbols())
        return symbol
    def get_by_ase(self):
        self.atoms=ase.io.read(self.__file_name__)
        self.bcell=self.atoms.get_reciprocal_cell()

