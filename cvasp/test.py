#!/usr/bin/python3
from plot_band import *
EIG,OCC,PRO,KPOINTS=plot_band.from_procar()
kpath,kpt_label,KPOINTS=plot_band.from_kpoints(path='.',read=False,save=True)
type1_Ele=[[["S"],'y','s'],[["Cu"],'r','>'],[["Sn"],'g','v'],[["Zn"],'b','<']]
figsize=False

EIG=plot_band.set_fermi(EIG,OCC)
plot_band.plot_band(kpath,kpt_label,EIG,PRO=PRO,type1_Ele=type1_Ele,type2_Ele_A=['Li','Cl'],figsize=figsize,plot_type=1)

