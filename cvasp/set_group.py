import re
import numpy as np
import plot_band
from plot_band import *
def set_group(PRO,SYMBOL,group):
    PRO=PRO[:,:,:,:,:10]
    promap={'s':[0],'py':[1],'pz':[2],'px':[3],'dxy':[4],'dyz':[5],'dz2':[6],'dx2-y2':[7],'dxz':[8],'p':[1,2,3],'d':[4,5,6,7,8]}
    PRO_GROUP=[]
    for igroup in group:
        tgroup=0
        #print(igroup)
        for iatom in igroup.split(','):
            if "_" in iatom:
                orbit=promap[iatom.split('_')[1]]
                whichatom=iatom.split('_')[0]
                if re.search('[A-Z]',whichatom):
                    whichatom=np.array([ i == whichatom for i in SYMBOL ]).nonzero()[0]
                else:
                    if re.search('-',whichatom):
                        whichatom=np.arange(int(whichatom.split('-')[0])-1,int(whichatom.split('-')[1]))
                    else:
                        whichatom=np.array([int(whichatom.split('-')[0])])
                tgroup=tgroup+PRO[:,:,:,whichatom][:,:,:,:,orbit].reshape((PRO.shape[0],PRO.shape[1],PRO.shape[2],len(whichatom),len(orbit))).sum((3,4))
            else:
                whichatom=iatom.split('_')[0]
                if re.search('[A-Z]',whichatom):
                    whichatom=np.array([ i == whichatom for i in SYMBOL ]).nonzero()[0]
                else:
                    if re.search('-',whichatom):
                        whichatom=np.arange(int(whichatom.split('-')[0])-1,int(whichatom.split('-')[1]))
                    else:
                        whichatom=np.array([int(whichatom.split('-')[0])-1])
                #print(whichatom)
                tgroup=tgroup+PRO[:,:,:,whichatom].reshape((PRO.shape[0],PRO.shape[1],PRO.shape[2],len(whichatom),PRO.shape[4])).sum((3,4))
        PRO_GROUP.append(tgroup)
    #print(len(PRO_GROUP))
    PRO_GROUP=np.array(PRO_GROUP)
    PRO_GROUP=PRO_GROUP/PRO.sum((3,4))
    print(PRO_GROUP.shape)
    return PRO_GROUP

##group=['N_pz','C_px,1_p,2','3']
#SYMBOL=plot_band.get_symbol(CONTCAR_file='./CONTCAR')
#PRO,EIG,OCC,TTOT,KPOINTS=plot_band.from_procar(ISPIN=1,read=False)
#PRO_GROUP=set_group(PRO,SYMBOL,group)
#print(PRO_GROUP.shape)
