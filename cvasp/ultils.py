import numpy as np
import re


def get_fermi(EIG,OCC,OCC_E=0.001):
    """
    this method used to get the fermi level through the occupid number, the occupid threshold are defined by OCC_E, default is 0.001.
    """
    FERMI=np.max(EIG[OCC>OCC_E])
    return FERMI

def set_group(project_data,grouptag,symbollist,norm=True):
    """
    this method used to  group the project data by given element|symbol|atomic number and/or orbital.
    project_data could be object:dos,pro from from_doscar.get_doscar() or from_procar.get_procar()
    grouptag should be writen as:
    symbollist shold be writen as string list.
    """
    if project_data.L_orbit ==  10:
        promap={'s':[0],'py':[1],'pz':[1],'px':[1],'dxy':[2],'dyz':[2],'dz2':[2],'dx2-y2':[2],'dxz':[2],'p':[1],'d':[2]}
        print('Lorbit=10')
    elif project_data.L_orbit > 10:
        promap={'s':[0],'py':[1],'pz':[2],'px':[3],'dxy':[4],'dyz':[5],'dz2':[6],'dx2-y2':[7],'dxz':[8],'p':[1,2,3],'d':[4,5,6,7,8]}
    project_=project_data.project
    pro_group=[]
    for igrouptag in grouptag:
        tot=0
        for iatom in igrouptag.split(','):
            whichorbit=[]
            if "_" in iatom:
                if bool(re.search('[xyz]',iatom.split('_')[1])) & (project_data.L_orbit ==  10):
                    print("Warnning: you are trying to project the band to lm orbit while using LORBIT==10")
                whichorbit=np.array(promap[iatom.split('_')[1]])
                whichatom=iatom.split('_')[0]
            else:
                whichatom=iatom
            if re.search('[A-Z]',whichatom): # element symbolliste e.g Ag,H,C ...
                whichatom=np.array([ i == whichatom for i in symbollist ]).nonzero()[0]
            else:
                if re.search('-',whichatom): # element number e.g 1-10,11-23 ...
                    whichatom=np.arange(int(whichatom.split('-')[0])-1,int(whichatom.split('-')[1]))
                else:
                    whichatom=np.array([int(whichatom)])-1 # element number e.g 1,2,3...
            if len(whichorbit) > 0 :
                if len(project_.shape) == 5:
                    tot += project_[:,:,:,whichatom,whichorbit].reshape((*project_.shape[0:3],whichatom.size,whichorbit.size)).sum((-1,-2))
                else: 
                    tot += project_[:,:,whichatom][:,:,:,whichorbit].reshape((*project_.shape[0:2],whichatom.size,whichorbit.size)).sum((-1,-2))
            else:
                if len(project_.shape) == 5:
                    tot += project_[:,:,:,whichatom].reshape((*project_.shape[0:3],whichatom.size,project_.shape[4])).sum((-1,-2))
                else:
                    tot += project_[:,:,whichatom].reshape((*project_.shape[0:2],whichatom.size,project_.shape[3])).sum((-1,-2))
        pro_group.append(tot)
    pro_group=np.array(pro_group,dtype=float)
    #if norm:
    #    print(pro_group.max())
    #    pro_group[:,project_.sum((-1,-2))!=0]=pro_group[:,project_.sum((-1,-2))!=0]/(project_.sum((-1,-2))[project_.sum((-1,-2))!=0])
    #    print(pro_group.max())
    return pro_group
