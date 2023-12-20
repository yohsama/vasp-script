#!/usr/bin/env python3
import ase.io
import ase
import sys
import numpy as np

pos=ase.Atoms()
if len(sys.argv)>1:
    pos=ase.io.read(sys.argv[1])
    prec=np.array(sys.argv[2],dtype='double')
else:
    pos=ase.io.read('POSCAR')
    prec=np.array(0.01)


rep=pos.cell.reciprocal() # get reciprocal vectors.
b1=np.array(rep[0])
b2=np.array(rep[1])
b3=np.array(rep[2])
cell=pos.get_cell()
a1=np.array(cell[0])
a2=np.array(cell[1])
a3=np.array(cell[2])
lattice_a,lattice_b,lattice_c,alpha,beta,gamma=pos.get_cell_lengths_and_angles()
alpha=alpha/180*np.pi
beta=beta/180*np.pi
gamma=gamma/180*np.pi

## all 14 kinds High-throughput electronic band structure.
## 10.1016/j.commatsci.2010.05.010
## High-throughput electronic band structure calculations: Challenges and tools
## Setyawaeta, Wahyu, Curtarolo, Stefano JCMS 2010

HighsymmK={}

def printkpath(HighsymmK,kpath):
    for i in kpath:
        if i=='|':
            print("---------")
        else:
            print(HighsymmK[i])
    return


##  A
##      A.1 Cubic (CUB, cP)
##          a1=(a, 0, 0)
##          a2=(0, a, 0)
##          a3=(0, 0, a)
##
a=cell.trace()/3
cell_stantard=np.diag([a,a,a])
if (np.abs(cell-cell_stantard)<prec).all():
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['M']      =[1/2,  1/2,    0]
    HighsymmK['R']      =[1/2,  1/2,  1/2]
    HighsymmK['X']      =[0,    1/2,    0]
    kpath=['Gamma','X','M','Gamma','X','|','M','R']
    print('CUB\n');printkpath(HighsymmK,kpath)
##      A.2 Face-centered cubic (FCC,cF)
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(0, a/2, a/2)               
##          a2=(0, a, 0)                  a2=(a/2, 0, a/2)
##          a3=(0, 0, a)                  a3=(a/2, a/2, 0)

a=(cell*np.array([[0,1,1],[1,0,1],[1,1,0]])).sum()/6*2
cell_stantard=np.array([[0, a/2, a/2],
                        [a/2, 0, a/2],
                        [a/2, a/2, 0]])

if (np.abs(cell-cell_stantard)<prec).all():
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['K']      =[3/8,  3/8,    3/4]
    HighsymmK['L']      =[1/2,  1/2,    1/2]
    HighsymmK['U']      =[5/8,  1/4,    5/8]
    HighsymmK['W']      =[1/2,  1/4,    3/4]
    HighsymmK['X']      =[1/2,    0,    1/2]
    kpath=['Gamma','X','W','K','Gamma','L','U','W','L','K','|','U','X']
    print('FCC\n');printkpath(HighsymmK,kpath)


##      A.3 Body-centered cubic (BCC, cI)
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(-a/2, a/2, a/2)               
##          a2=(0, a, 0)                  a2=(a/2, -a/2, a/2)
##          a3=(0, 0, a)                  a3=(a/2, a/2, -a/2)
##

a=(cell*np.array([[-1,1,1],[1,-1,1],[1,1,-1]])).sum()/6*2
cell_stantard=np.array([[-a/2, a/2, a/2],
                        [a/2, -a/2, a/2],
                        [a/2, a/2, -a/2]])

if (np.abs(cell-cell_stantard)<prec).all():
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['H']      =[1/2,  -1/2,    1/2]
    HighsymmK['P']      =[1/4,  1/4,    1/4]
    HighsymmK['N']      =[0,  0,    1/2]
    kpath=['Gamma','H','N','Gamma','P','H','|','P','N']
    print('BCC\n');printkpath(HighsymmK,kpath)


##      A.4 Tetragonal (TET, tP) Lattice
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(-a/2, a/2, a/2)               
##          a2=(0, a, 0)                  a2=(a/2, -a/2, a/2)
##          a3=(0, 0, a)                  a3=(a/2, a/2, -a/2)
##

a=(cell*np.array([[1,0,0],[0,1,0],[0,0,0]])).sum()/2
c=cell[2,2]
cell_stantard=np.array([[a, 0, 0],
                        [0, a, 0],
                        [0, 0, c]])

if (np.abs(cell-cell_stantard)<prec).all():
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['A']      =[1/2,  1/2,    1/2]
    HighsymmK['M']      =[1/2,  1/2,    0]
    HighsymmK['R']      =[0,    1/2,    1/2]
    HighsymmK['X']      =[0,    1/2,    0]
    HighsymmK['Z']      =[0,    0,    1/2]
    kpath=['Gamma','X','M','Gamma','Z','R','A','Z','|','X','R','|','M','A']
    print('TET\n');printkpath(HighsymmK,kpath)






##      A.5. Body-centered tetragonal (BCT, tI) 
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(-a/2, a/2, c/2)               
##          a2=(0, a, 0)                  a2=(a/2, -a/2, c/2)
##          a3=(0, 0, c)                  a3=(a/2, a/2, -c/2)


a=(cell*np.array([[-1,1,0],[1,-1,0],[1,1,0]])).sum()/6*2
c=(cell*np.array([[0,0,1],[0,0,1],[0,0,-1]])).sum()/3*2
cell_stantard=np.array([[-a/2, a/2, c/2],
                        [a/2, -a/2, c/2],
                        [a/2, a/2, -c/2]])

if (np.abs(cell-cell_stantard)<prec).all():
    # there are 2 varialtions in BCT
    if c<a:                #BCT1
        eta=(1+c**2/a**2)/4
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['M']      =[-1/2,  1/2,    1/2]
        HighsymmK['N']      =[0,    1/2,    0]
        HighsymmK['P']      =[1/4,    1/4,    1/4]
        HighsymmK['X']      =[0,    0,    1/2]
        HighsymmK['Z']      =[eta,    eta,      -eta]
        HighsymmK['Z_1']      =[-eta,    1-eta,    eta]
        kpath=['Gamma','X','M','Gamma','Z','P','N','Z_1','M','|','X','P']
        print('BCT1\n');printkpath(HighsymmK,kpath)

    else:                #BCT2
        eta=(1+c**2/a**2)/4
        ksi=a**2/(2*c**2)
        HighsymmK['Gamma']  =[0 ,       0,      0]
        HighsymmK['N']      =[0,        1/2,    0]
        HighsymmK['P']      =[1/4,      1/4,    1/4]
        HighsymmK['Sigma']  =[-eta,       eta,      eta]
        HighsymmK['Sigma_1']=[eta,       1-eta,      -eta]
        HighsymmK['X']      =[0,        0,      1/2]
        HighsymmK['Y']      =[-ksi,      ksi,      1/2]
        HighsymmK['Y_1']      =[1/2,      1/2,      -ksi]
        HighsymmK['Z']      =[1/2,        1/2,      -1/2]
        kpath=['Gamma','X','Y','Sigma','Gamma','Z','Sigma_1','N','P','Y_1','Z','|','X','P']
        print('BCT2\n');printkpath(HighsymmK,kpath)




##      A.6. Orthorhombic (ORC, oP)   
##          a1=(a, 0, 0)                              
##          a2=(0, b, 0)               
##          a3=(0, 0, c)               



a=cell[0,0]
b=cell[1,1]
c=cell[2,2]
cell_stantard=np.array([[a,0,0],
                        [0,b,0],
                        [0,0,c]])
if not a<b<c:
    print("Orthorhombic request a<b<c")
    


if (np.abs(cell-cell_stantard)<prec).all():
    print((np.abs(cell-cell_stantard)<prec))
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['R']      =[1/2,    1/2,    1/2]
    HighsymmK['S']      =[1/2,  1/2,    0]
    HighsymmK['U']      =[1/2,  0,    1/2]
    HighsymmK['T']      =[0,    1/2,    1/2]
    HighsymmK['X']      =[1/2,    0,    0]
    HighsymmK['Y']      =[0,    1/2,    0]
    HighsymmK['Z']      =[0,    0,    1/2]
    kpath=['Gamma','X','S','Y','Gamma','Z','U','R','T','Z','|','Y','T','|','U','X','|','S','R']
    print('ORC\n');printkpath(HighsymmK,kpath)


##      A.7. Face-centered orthorhombic (ORCF, oF) Ordering
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(0, b/2, c/2)               
##          a2=(0, a, 0)                  a2=(a/2, 0, c/2)
##          a3=(0, 0, c)                  a3=(a/2, b/2, 0)

a=(cell[1,0]+cell[2,0])*2
b=(cell[0,1]+cell[2,1])*2
c=(cell[0,2]+cell[1,2])*2
cell_stantard=np.array([[0, b/2, c/2],
                        [a/2, 0, c/2],
                        [a/2, b/2, 0]])

if (np.abs(cell-cell_stantard)<prec).all():
    # there are 2 varialtions in ORCF
    if (1/a**2>=(1/b**2+1/c**2)):       #ORCF1,OCF3 when reach "="
        eta=(1+a**2/b**2+a**2/c**2)/4
        ksi=(1+a**2/b**2-a**2/c**2)/4
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['A']      =[1/2,  1/2+ksi,    ksi]
        HighsymmK['A_1']    =[1/2,    1/2-ksi,    1-ksi]
        HighsymmK['L']      =[1/2,    1/2,    1/2]
        HighsymmK['T']      =[1,    1/2,    1/2]
        HighsymmK['X']      =[0,    eta,      eta]
        HighsymmK['X_1']    =[1,    1-eta,    1-eta]
        HighsymmK['Y']      =[1/2,    0,      1/2]
        HighsymmK['Z']      =[1/2,    1/2,    0]
        kpath=['Gamma','Y','T','Z','Gamma','X','A_1','Y','|','T','X_1','|','X','A','Z','|','L','Gamma']
        print('ORCF\n');printkpath(HighsymmK,kpath)

    else:               #OCF2
        eta=(1+a**2/b**2-a**2/c**2)/4
        delta=(1+b**2/a**2-b**2/c**2)/4
        phi=(1+c**2/b**2-c**2/a**2)/4
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['C']      =[1/2,  1/2-eta,    1-eta]
        HighsymmK['C_1']    =[1/2,    1/2+eta,    eta]
        HighsymmK['D']      =[1/2-delta,  1/2,    1-delta]
        HighsymmK['D_1']    =[1/2+delta,    1/2,    delta]
        HighsymmK['L']      =[1/2,    1/2,    1/2]
        HighsymmK['H']      =[1-phi,    1/2-phi,    1/2]
        HighsymmK['H_1']    =[phi,    1/2+phi,      1/2]
        HighsymmK['X']      =[0,    1/2,    1/2]
        HighsymmK['Y']      =[1/2,    0,      1/2]
        HighsymmK['Z']      =[1/2,    1/2,    0]
        kpath=['Gamma','Y','C','D','X','Gamma','Z','D_1','H','C','|','C_1','Z','|','X','AH_1','|','H','Y','|','L','Gamma']
        print('ORCF\n');printkpath(HighsymmK,kpath)




##      A.8. Body-centered orthorhombic (ORCI, oI) Ordering
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(-a/2, b/2, c/2)               
##          a2=(0, a, 0)                  a2=(a/2, -b/2, c/2)
##          a3=(0, 0, c)                  a3=(a/2, b/2, -c/2)


a=(-cell[0,0]+cell[1,0]+cell[2,0])/3*2
b=( cell[0,1]-cell[0,1]+cell[2,1])/3*2
c=( cell[0,2]+cell[1,2]-cell[2,2])/3*2
cell_stantard=np.array([[-a/2,b/2,c/2],
                        [a/2,-b/2,c/2],
                        [a/2,b/2,-c/2]])

if not a<b<c:
    print("ORCI request a<b<c")
    



if (np.abs(cell-cell_stantard)<prec).all():
    ksi=(1+a**2/c**2)/4
    eta=(1+b**2/c**2)/4
    delta=(b**2-a**2)/(4*c**2)
    miu=(a**2+b**2)/(4*c**2)
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['L']      =[-miu ,   miu,      1/2-delta]
    HighsymmK['L_1']    =[ miu ,  -miu,      1/2+delta]
    HighsymmK['L_2']    =[1/2-delta ,  1/2+delta, -miu]
    HighsymmK['R']      =[0,    1/2,    0]
    HighsymmK['S']      =[1/2,    0,    0]
    HighsymmK['T']      =[0,    0,    1/2]
    HighsymmK['W']      =[1/4,    1/4,    1/4]
    HighsymmK['X']      =[-ksi,    ksi,    ksi]
    HighsymmK['X_1']    =[ksi,    1-ksi,   -ksi]
    HighsymmK['Y']      =[eta,    -eta,    eta]
    HighsymmK['Y_1']    =[1-eta,    eta,   -eta]
    HighsymmK['Z']      =[1/2,    1/2,    -1/2]
    kpath=['Gamma','X','L','T','W','R','X_1','Z','Gamma','Y','S','W','|','L_1','Y','|','Y_1','Z']
    print('ORCI\n');printkpath(HighsymmK,kpath)


##      A.9. C-centered orthorhombic (ORCC, oS) Ordering
##          Conventional lattice        Primitive lattice
##          a1=(a, 0, 0)                  a1=(a/2, -b/2, 0)               
##          a2=(0, a, 0)                  a2=(a/2, b/2, 0)
##          a3=(0, 0, c)                  a3=(0 , 0, c)


a=( cell[0,0]-cell[1,0])/2*2
b=( cell[0,1]+cell[0,1])/2*2
c=( cell[2,2])
cell_stantard=np.array([[a/2,-b/2,0],
                        [a/2,b/2,0],
                        [0,0,c]])


if (np.abs(cell-cell_stantard)<prec).all():
    ksi=(1+a**2/b**2)/4
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['A']      =[ksi,  ksi,    1/2]
    HighsymmK['A_1']    =[-ksi, 1-ksi,   1/2]
    HighsymmK['R']      =[0,    1/2,    1/2]
    HighsymmK['S']      =[0,    1/2,    0]
    HighsymmK['T']      =[-1/2,    1/2,    1/2]
    HighsymmK['X']      =[ksi,    ksi,    0]
    HighsymmK['X_1']    =[-ksi,    1-ksi,   0]
    HighsymmK['Y']      =[-1/2,    1/2,    0]
    HighsymmK['Z']      =[0,    0,    1/2]
    kpath=['Gamma','X','S','R','A','Z','Gamma','Y','X_1','A_1','T','Y','|','Z','T']
    print('ORCC\n');printkpath(HighsymmK,kpath)



##      A.10. Hexagonal (HEX, hP) Lattice
##          a1=(a/2, -(a*sqrt(3)/2, 0)               
##          a2=(a/2, (a*sqrt(3))/2, 0)
##          a3=(0 , 0, c)


a=  (( cell[0,0]-cell[1,0]/np.sqrt(3))*2 +
    ( cell[0,1]+cell[0,1]/np.sqrt(3))*2)/4
c=  ( cell[2,2])
cell_stantard=np.array([[a/2,-a*np.sqrt(3)/2,0],
                        [a/2,a*np.sqrt(3)/2,0],
                        [0,0,c]])


if (np.abs(cell-cell_stantard)<prec).all():
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['A']      =[0,    0,      1/2]
    HighsymmK['H']      =[1/3,  1/3,    1/2]
    HighsymmK['K']      =[1/3,  1/3,    0]
    HighsymmK['L']      =[1/2,    0,    1/2]
    HighsymmK['M']      =[1/2,    0,    0]
    kpath=['Gamma','M','K','Gamma','A','L','H','A','|','L','M','|','K','H']


##      A.11. Rhombohedral (RHL, hR) Lattice
##          a1=(a*cos(alpha/2), -a*sin(alpha/2), 0)               
##          a2=(a*cos(alpha/2), a*sin(alpha/2), 0)
##          a3=(a*cos(alpha)/cos(alpha/2), 0, a*sqrt(1-cos^2(alpha)/cos^2(alpha/2)) )

a=np.mean((lattice_a,lattice_b,lattice_c))

#a = (np.sqrt(cell[0,0]**2+c[0,1]**2)+np.sqrt(cell[1,0]**2+c[1,1]**2))/2

cell_stantard=np.array([[a*np.cos(alpha/2),-a*np.sin(alpha/2),0],
                        [a*np.cos(alpha/2),-a*np.sin(alpha/2),0],
                        [a*np.cos(alpha)/np.cos(alpha/2),0,a*np.sqrt(1-np.cos(alpha)**2/np.cos(alpha/2)**2)]])


if (np.abs(cell-cell_stantard)<prec).all():
    # there are 2 varialtions in RHL
    if alpha<np.pi/2:    #RHL1
        eta=(1+4*np.cos(alpha)/(2+4*np.cos(alpha)))
        niu=3/4-eta/2
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['B']      =[eta,    1/2,   1-eta]
        HighsymmK['B_1']    =[1/2,    1-eta,  eta-1]
        HighsymmK['F']      =[1/2,  1/2,    0]
        HighsymmK['L']      =[1/2,    0,    0]
        HighsymmK['L_1']      =[0,    0,    -1/2]
        HighsymmK['P']      =[0,    0,    -1/2]
        HighsymmK['P_1']      =[0,    0,    -1/2]
        HighsymmK['P_2']      =[0,    0,    -1/2]
        HighsymmK['Q']      =[1-niu,    niu,    ]
        HighsymmK['X']      =[niu,    0,    -niu]
        HighsymmK['Z']      =[1/2,    1/2,    1/2]
        kpath=['Gamma','L','B_1','|','B','Z','Gamma','X','Q','F','P_1','Z','|','L','P']
        print('RHL1\n');printkpath(HighsymmK,kpath)

    else:    #RHL2
        eta=1/2*(np.tan(alpha/2)**2)
        niu=3/4-eta/2
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['F']      =[1/2,    1/2,   0]
        HighsymmK['L']    =[1/2,    0,  0]
        HighsymmK['P']      =[1-niu,  -niu,    1-niu]
        HighsymmK['P_1']      =[niu,    niu-1,    niu-1]
        HighsymmK['Q']      =[eta,    eta,    eta]
        HighsymmK['Q_1']      =[1-eta,    -eta,    -eta]
        HighsymmK['Z']      =[1/2,    -1/2,    -1/2]
        kpath=['Gamma','L','P','Z','Q','Gamma','F','P_1','Q_1','L','Z']
        print('RHL2\n');printkpath(HighsymmK,kpath)



##      A.12. Monoclinic (MCL, mP)
##          a1=(a*cos(alpha/2), -a*sin(alpha/2), 0)               
##          a2=(a*cos(alpha/2), a*sin(alpha/2), 0)
##          a3=(a*cos(alpha)/cos(alpha/2), 0, a*sqrt(1-cos^2(alpha)/cos^2(alpha/2)) )

a=cell[0,0]
b=cell[1,1]
c=lattice_c

if not (a<=c)&(b<=c)&(alpha<(np.pi/2)):
    print("MCL request a,b<=c & alpha<90")

#a = (np.sqrt(cell[0,0]**2+c[0,1]**2)+np.sqrt(cell[1,0]**2+c[1,1]**2))/2

cell_stantard=np.array([[a,0,0],
                        [0,b,0],
                        [0,c*np.cos(alpha),c*np.sin(alpha)]])

if (np.abs(cell-cell_stantard)<prec).all():
    eta=(1-b*np.cos(alpha)/c/(2*np.sin(alpha)**2))
    niu=1/2-eta*c*np.cos(alpha)/b
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['A']      =[1/2,    1/2,   0]
    HighsymmK['C']      =[0,    1/2,   1/2]
    HighsymmK['D']      =[1/2,    0,   1/2]
    HighsymmK['D_1']    =[1/2,    0,    -1/2]
    HighsymmK['E']      =[1/2,    1/2,    1/2]
    HighsymmK['H']      =[0,    eta,     1-niu]
    HighsymmK['H_1']      =[0,    1-eta,   niu]
    HighsymmK['H_2']      =[0,    eta,     niu]
    HighsymmK['M']      =[1/2,    eta,    1-niu]
    HighsymmK['M_1']      =[1/2,    1-eta,   niu ]
    HighsymmK['M_2']      =[1/2,    eta,    -niu]
    HighsymmK['X']      =[0,    1/2,    0]
    HighsymmK['Y']      =[0,    0,   1/2]
    HighsymmK['Y_1']      =[0,    0,    -1/2]
    HighsymmK['Z']      =[1/2,    0,    0]
    kpath=['Gamma','Y','H','C','E','M_1','A','X','H_1','|','M','D','Z','|','Y','D']
    print('MCL\n');printkpath(HighsymmK,kpath)




##      A.13. C-centered monoclinic (MCLC, mS) 
##          a1=(a*cos(alpha/2), -a*sin(alpha/2), 0)               
##          a2=(a*cos(alpha/2), a*sin(alpha/2), 0)
##          a3=(a*cos(alpha)/cos(alpha/2), 0, a*sqrt(1-cos^2(alpha)/cos^2(alpha/2)) )

a=(cell[0,0]-cell[1,0])
b=(cell[0,1]+cell[1,1])
c=np.sqrt(cell[2,1]**2+cell[2,2]**2)
calpha=np.arcsin(cell[2,2]/c)
if not (a<=c)&(b<=c)&(alpha<(np.pi/2)):
    print(lattice_a,lattice_b,lattice_c,alpha)
    print(a,b,c)
    print("MCL request a,b<=c & alpha<90")

#a = (np.sqrt(cell[0,0]**2+c[0,1]**2)+np.sqrt(cell[1,0]**2+c[1,1]**2))/2

cell_stantard=np.array([[a/2,b/2,0],
                        [-a/2,b/2,0],
                        [0,c*np.cos(calpha),c*np.sin(calpha)]])

print(cell_stantard)
if (np.abs(cell-cell_stantard)<prec).all():
    # there are 2 varialtions in ORCF
    k_gamma=np.arccos(b1.dot(b2)/np.linalg.norm(b1)/np.linalg.norm(b2))
    if (k_gamma-np.pi/2)>-prec:           #MCLC1 &MCLC2
        ksi=(2-b*np.cos(alpha)/c)/(4*np.sin(alpha)**2)
        eta=1/2+(2*ksi*c*np.cos(alpha)/b)
        psi=3/4-alpha**2/(4*b**2*np.sin(alpha))
        phi=psi+(3/4-psi)*b*np.cos(alpha)/c
        HighsymmK['Gamma']  =[0 ,   0,      0]
        HighsymmK['N']      =[1/2,    0,   0]
        HighsymmK['N_1']      =[0,    1/2,   0]
        HighsymmK['F']      =[1-ksi,    1-ksi,   1-eta]
        HighsymmK['F_1']    =[ksi,    ksi,    eta]
        HighsymmK['F_2']      =[-ksi,    -ksi,    1-eta]
        HighsymmK['F_3']      =[1-ksi,    -ksi,     1-eta]
        HighsymmK['I']      =[phi,    1-phi,  1/2]
        HighsymmK['I_1']      =[1-phi,    phi-1,     1/2]
        HighsymmK['L']      =[1/2,    1/2,    1/2]
        HighsymmK['M']      =[1/2,    0,   1/2]
        HighsymmK['X']      =[1-psi,    psi-1,    0]
        HighsymmK['X_1']      =[psi,    1-psi,    0]
        HighsymmK['X_2']      =[psi-1,    -psi,   0]
        HighsymmK['Y']      =[1/2,    1/2,    0]
        HighsymmK['Y_1']      =[-1/2,    -1/2,    0]
        HighsymmK['Z']      =[0,    0,    1/2]
        if (k_gamma-np.pi/2) > prec:          #MCLC1 
            kpath=['Gamma','Y','F','L','I','|','I_1','Z','F_1','|','Y','X_1','|','X','Gamma','N','|','M','Gamma']
            print('MCLC1\n');printkpath(HighsymmK,kpath)
        else:                               #MCLC2
            kpath=['Gamma','Y','F','L','I','|','I_1','Z','F_1','|','N','Gamma','M']
            print('MCLC2\n');printkpath(HighsymmK,kpath)

    else:            #MCLC3-5:
        if (b*np.cos(alpha)/c+b**2*np.sin(alpha)**2/a**2 - 1) < prec:  #MCLC3 & MCLC4  
            miu=(1+b**2/a**2)/4
            delta=b*c*np.cos(alpha)/(2*a**2)
            ksi=miu-1/4+(1-b*np.cos(alpha)/c)/(4*np.sin(alpha)**2)
            eta=1/2+(2*ksi*c*np.cos(alpha)/b)
            phi=1+ksi-2*miu
            psi=eta-2*delta
            HighsymmK['Gamma']  =[0 ,   0,      0]
            HighsymmK['F']      =[1-phi,    1-phi,   1-psi]
            HighsymmK['F_1']    =[phi,    phi-1,    psi]
            HighsymmK['F_2']      =[1-phi,    -phi,    1-psi]
            HighsymmK['H']      =[ksi,    ksi,     eta]
            HighsymmK['H_1']      =[1-ksi,    -ksi,  1-eta]
            HighsymmK['H_2']      =[ -ksi,    -ksi,  1-eta]
            HighsymmK['I']      =[1/2,    -1/2,    1/2]
            HighsymmK['M']      =[1/2,    0,   1/2]
            HighsymmK['N']      =[1/2,    0,    0]
            HighsymmK['N_1']      =[0,    -1/2,    0]
            HighsymmK['X']      =[1/2,    -1/2,   0]
            HighsymmK['Y']      =[miu,    miu,    delta]
            HighsymmK['Y_1']      =[1-miu,    -miu,    -delta]
            HighsymmK['Y_2']      =[-miu,    -miu,    -delta]
            HighsymmK['Y_3']      =[miu,    miu-1,    delta]
            HighsymmK['Z']      =[0,    0,    1/2]
            if (b*np.cos(alpha)/c+b**2*np.sin(alpha)**2/a**2 - 1) < -prec: #MCLC3
                kpath=['Gamma','Y','F','H','Z','I','F_1','|','H_1','Y_1','X','Gamma','N','|','M','Gamma']
                print('MCLC3\n');printkpath(HighsymmK,kpath)
            if (b*np.cos(alpha)/c+b**2*np.sin(alpha)**2/a**2 - 1) < -prec: #MCLC3
                kpath=['Gamma','Y','F','H','Z','I','|','H_1','Y_1','X','Gamma','N','|','M','Gamma']
                print('MCLC4\n');printkpath(HighsymmK,kpath)
        else:  #MCLC5
            ksi=(b**2/a**2+(1-b*np.cos(alpha)/c)/np.sin(alpha)**2)/4
            eta=1/2+(2*ksi*c*np.cos(alpha)/b)
            miu=eta/2+b**2/(4*a**2)-b*c*np.cos(alpha)/(2*a**2)
            niu=2*miu-ksi
            omega=(4*niu-1-b**2*np.sin(alpha)**2/a**2)*c/(2*b*np.cos(alpha))
            delta=ksi*c*np.cos(alpha)/b+omega/2-1/4
            pho=1-ksi*a**2/b**2
            HighsymmK['Gamma']  =[0 ,   0,      0]
            HighsymmK['F']      =[niu,  niu,    omega]
            HighsymmK['F_1']    =[1-niu,1-niu,1-omega]
            HighsymmK['F_2']      =[niu,niu-1,omega]
            HighsymmK['H']      =[ksi,ksi,eta]
            HighsymmK['H_1']      =[1-ksi,    -ksi,  1-eta]
            HighsymmK['H_2']      =[ -ksi,    -ksi,  1-eta]
            HighsymmK['I']      =[pho,1-pho,1/2]
            HighsymmK['I_1']      =[1-pho,pho-1,1/2]
            HighsymmK['L']      =[1/2,    1/2,    1/2]
            HighsymmK['M']      =[1/2,    0,    1/2]
            HighsymmK['N']      =[1/2,    0,    0]
            HighsymmK['N_1']      =[0,    -1/2,    0]
            HighsymmK['X']      =[1/2,    -1/2,   0]
            HighsymmK['Y']      =[miu,    miu,    delta]
            HighsymmK['Y_1']      =[1-miu,    -miu,    -delta]
            HighsymmK['Y_2']      =[-miu,    -miu,    -delta]
            HighsymmK['Y_3']      =[miu,    miu-1,    delta]
            HighsymmK['Z']      =[0,    0,    1/2]
            kpath=['Gamma','Y','F','L','I','|','I_1','Z','H','F_1','|','H_1','Y_1','X','Gamma','N','|','M','Gamma']
            print('MCLC5\n');printkpath(HighsymmK,kpath)


##      A.14. Triclinic (TRI, aP) Lattice
##          a1=            
##          a2=
##          a3=

k_gamma=np.arccos(b1.dot(b2)/np.linalg.norm(b1)/np.linalg.norm(b2))
k_alpha=np.arccos(b2.dot(b3)/np.linalg.norm(b3)/np.linalg.norm(b2))
k_beta=np.arccos(b3.dot(b1)/np.linalg.norm(b1)/np.linalg.norm(b3))
print('kg=%f,kb=%f,ka=%f' % (k_gamma/np.pi*180,k_beta/np.pi*180,k_alpha/np.pi*180))
if not (((np.pi/2-prec)<=k_gamma<k_alpha)&((np.pi/2-prec)<=k_gamma<k_beta)):
    print("TRI request \n")
    print(" TRI_1a k_alpha > k_gamma, k_beta > k_gamma, k_gamma >= 90, \n")

if not (((np.pi/2+prec)>=k_gamma>k_alpha)&((np.pi/2+prec)>=k_gamma>k_beta)):
    print("TRI request \n")
    print(" TRI_1b k_alpha < k_gamma, k_beta < k_gamma, k_gamma <= 90, \n")

if (k_gamma-np.pi/2) > -prec:    #TRI_1a &   TRI_2a
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['L']      =[1/2,    1/2,    0]
    HighsymmK['M']      =[0,  1/2,    1/2]
    HighsymmK['N']      =[1/2,  0,    1/2]
    HighsymmK['R']      =[1/2,    1/2,    1/2]
    HighsymmK['X']      =[1/2,    0,    0]
    HighsymmK['Y']      =[0,    1/2,    0]
    HighsymmK['Z']      =[0,    0,    1/2]
    kpath=['X','Gamma','Y','|','L','Gamma','Z','|','N','Gamma','M','|','R','Gamma']
    print('TRI1\n');printkpath(HighsymmK,kpath)
elif (k_gamma-np.pi/2) < prec:    #TRI_1b &   TRI_2b
    HighsymmK['Gamma']  =[0 ,   0,      0]
    HighsymmK['L']      =[1/2,   -1/2,    0]
    HighsymmK['M']      =[0,  0,    1/2]
    HighsymmK['N']      =[-1/2,  -1/2,    1/2]
    HighsymmK['R']      =[0,    -1/2,    1/2]
    HighsymmK['X']      =[0,    -1/2,    0]
    HighsymmK['Y']      =[1/2,    0,    0]
    HighsymmK['Z']      =[-1/2,    0,    1/2]
    kpath=['X','Gamma','Y','|','L','Gamma','Z','|','N','Gamma','M','|','R','Gamma']
    print('TRI2\n');printkpath(HighsymmK,kpath)

