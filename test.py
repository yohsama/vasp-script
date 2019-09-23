#!/usr/bin/python3
from plot_band import *
import sys
ISPIN=int(input("ISPIN\n"))
EIG,OCC,PRO,KPOINTS=plot_band.from_procar(ISPIN=ISPIN,read=True)
kpath,kpt_label,KPOINTS=plot_band.from_kpoints(path='.',read=False,save=True)
type=int(input("fat_band/1,colormap_band/2\n"))
type1_Ele=[]
type2_Ele_A=[]
type2_Ele_B=""
if type==1:
    while True:
        tmp=(input("Element,split by \",\"\n"))
        Element=tmp.split(",")
        if Element==[]:
            exit()
        else:
            color=(input("color\n"))
            symbol=input("symbol\n")
            type1_Ele.append([Element,color,symbol])
        tmp=input("Add New Element? ENTER for EXIT, any words for CONTINUE\n")
        if tmp =="":
            break
if type==2:
    tmp=(input("Element,split by \",\"\n"))
    type2_Ele_A=tmp.split(",")
    tmp=input("Part B ? ENTER for EXIT, any words for CONTINUE\n")
    if not tmp =="":
        tmp=(input("Element,split by \",\"\n"))
        type2_Ele_B=tmp.split(",")


Elim=(int(input("Energy ymin\n")),int(input("Energy ymax\n")))
figsize=False
EIG=plot_band.set_fermi(EIG,OCC)
plot_band.plot_band(kpath,kpt_label,EIG,PRO=PRO,Elim=Elim,type1_Ele=type1_Ele,type2_Ele_A=type2_Ele_A,type2_Ele_B=type2_Ele_B,figsize=figsize,plot_type=type)

