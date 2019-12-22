#!/usr/bin/env python3
import numpy as np
import struct
import scipy
import scipy.integrate
import time
import sys
import ase.io
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from multiprocessing import Pool
path='./'
def get_wave(path,selectband):
    t1=time.time()
    SYS=ase.io.read(path+'/CONTCAR',format='vasp')
    wavefile=open(path+'/WAVECAR','rb')
    Recl=int(struct.unpack('d',wavefile.read(8))[0]/8)
    Ispin=int(struct.unpack('d',wavefile.read(8))[0])
    RFLAG=int(struct.unpack('d',wavefile.read(8))[0])
    wavefile.seek((Recl-3)*8,1)
    N_wave_k_vector=int(struct.unpack('d',wavefile.read(8))[0])
    N_Band=int(struct.unpack('d',wavefile.read(8))[0])
    Encut=(struct.unpack('d',wavefile.read(8))[0])
    cell=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cell[i,j]=struct.unpack('d',wavefile.read(8))[0]
    c=0.26246582250210965422
    b=np.zeros((3,3))
    Vcell=np.dot(np.cross(cell[0,:],cell[1,:]),cell[2,:])
    b[0,:]=2*np.pi*np.cross(cell[1,:],cell[2,:])/Vcell
    b[1,:]=2*np.pi*np.cross(cell[2,:],cell[0,:])/Vcell
    b[2,:]=2*np.pi*np.cross(cell[0,:],cell[1,:])/Vcell
    bmag=np.sqrt(np.sum(b**2,1))
    phi12=np.arccos(np.dot(b[0,:],b[1,:])/bmag[0]*bmag[1])
    vtmp=np.cross(b[0,:],b[1,:])
    vmag=np.linalg.norm(np.cross(b[0,:],b[1,:]))
    sinphi123=np.dot(b[2,:],vtmp/(vmag*bmag[2]))
    nb1maxA=np.floor(np.sqrt(Encut*c)/(bmag[0]*abs(np.sin(phi12)))+1)
    nb2maxA=np.floor(np.sqrt(Encut*c)/(bmag[1]*abs(np.sin(phi12)))+1)
    nb3maxA=np.floor(np.sqrt(Encut*c)/(bmag[2]*abs(sinphi123))+1)
    npmaxA=np.round(4*np.pi*nb1maxA*nb2maxA*nb3maxA/3)
    phi13=np.arccos(np.dot(b[0,:],b[2,:])/bmag[0]*bmag[2])
    vtmp=np.cross(b[0,:],b[2,:])
    vmag=np.linalg.norm(vtmp)
    sinphi123=np.dot(b[1,:],vtmp/(vmag*bmag[1]))
    nb1maxB=np.floor(np.sqrt(Encut*c)/bmag[0]*abs(np.sin(phi13))+1)
    nb2maxB=np.floor(np.sqrt(Encut*c)/bmag[1]*abs(sinphi123)+1)
    nb3maxB=np.floor(np.sqrt(Encut*c)/bmag[2]*abs(np.sin(phi13))+1)
    npmaxB=np.round(4*np.pi*nb1maxB*nb2maxB*nb3maxB/3)
    phi23=np.arccos(np.dot(b[1,:],b[2,:])/bmag[1]*bmag[2])
    vtmp=np.cross(b[1,:],b[2,:])
    vmag=np.linalg.norm(vtmp)
    sinphi123=np.dot(b[0,:],vtmp/(vmag*bmag[1]))
    nb1maxC=np.floor(np.sqrt(Encut*c)/bmag[0]*abs(sinphi123)+1)
    nb2maxC=np.floor(np.sqrt(Encut*c)/bmag[1]*abs(np.sin(phi23))+1)
    nb3maxC=np.floor(np.sqrt(Encut*c)/bmag[2]*abs(np.sin(phi23))+1)
    npmaxC=np.round(4*np.pi*nb1maxC*nb2maxC*nb3maxC/3)
    nb1max=int(np.max([nb1maxA,nb1maxB,nb1maxC]))
    nb2max=int(np.max([nb2maxA,nb2maxB,nb2maxC]))
    nb3max=int(np.max([nb3maxA,nb3maxB,nb3maxC]))
    npmax=int(np.max([npmaxA,npmaxB,npmaxC]))
    sumkg=np.zeros(3)
    igall=np.zeros((Ispin,N_wave_k_vector,Recl,3),dtype='double')
    coeff=np.zeros((Ispin,N_wave_k_vector,N_Band,Recl),dtype='complex')
    #coeff=np.zeros((Ispin,N_wave_k_vector,N_Band,Recl),dtype='double')
    Engvaule=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
    #eig=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
    eig=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
    occ=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
    wave_k_vector=np.zeros((Ispin,N_wave_k_vector,3),dtype='double')
    wavefile.seek(8*(Recl-12),1)
    print(time.time()-t1,"pre-set done, start read coeff")
    for ispin in range(Ispin):
        for ikpt in range(N_wave_k_vector):
            N_plane=int(struct.unpack('d',wavefile.read(8))[0])
            for i in range(3):
                wave_k_vector[ispin,ikpt,i]=struct.unpack('d',wavefile.read(8))[0]
            t1=time.time()
            eig_tmp=np.array(struct.unpack('d'*3*N_Band,wavefile.read(8*3*N_Band)));
            eig[ispin,ikpt]=np.sign(eig_tmp[0::3])*np.abs(eig_tmp[0::3]+eig_tmp[1::3]*1j)
            occ[ispin,ikpt]=eig_tmp[2::3]
            #for iband in range(N_Band):
            #    #eig[ispin,ikpt,iband]=struct.unpack('d',wavefile.read(8))[0]
            #    eig[ispin,ikpt,iband]=struct.unpack('d',wavefile.read(8))[0]+struct.unpack('d',wavefile.read(8))[0]*1j
            #    occ[ispin,ikpt,iband]=struct.unpack('d',wavefile.read(8))[0]
            print(time.time()-t1,"eig, occ done")
            ncnt=-1
            t1=time.time()
            nb1=list(range(0,nb1max))+list(range(-nb1max-1,0))
            nb2=list(range(0,nb2max))+list(range(-nb2max-1,0))
            nb3=list(range(0,nb3max))+list(range(-nb3max-1,0))
            igtmp=np.array(np.meshgrid(nb2,nb3,nb1)).reshape(3,(2*nb3max+1)*(2*nb2max+1)*(2*nb1max+1)).T
            igtmp=igtmp[:,np.array([2,0,1])]
            sumkg=np.matrix(igtmp+wave_k_vector[ispin,ikpt])*np.matrix(b)
            gtot=np.linalg.norm(sumkg,axis=1)
            etot=(gtot**2)/c
            igtmp=igtmp[etot<=Encut]
            igtmp=igtmp.dot(b)
            np.savetxt('igtmp',igtmp)
            igall[ispin,ikpt,0:igtmp.shape[0],:]=igtmp
           # for ig3p in list(range(0,nb3max))+list(range(-nb3max-1,0)):
           #     for ig2p in list(range(0,nb2max))+list(range(-nb2max-1,0)):
           #         for ig1p in list(range(0,nb1max))+list(range(-nb1max-1,0)):
           #             for j in range(3):
           #                 sumkg[j]=(wave_k_vector[ispin,ikpt,0]+ig1p)*b[0,j]+(wave_k_vector[ispin,ikpt,1]+ig2p)*b[1,j]+(wave_k_vector[ispin,ikpt,2]+ig3p)*b[2,j]
           #             gtot=np.linalg.norm(sumkg)
           #             etot=gtot**2/c
           #             if (etot<=Encut):
           #                 ncnt=ncnt+1
           #                 igall[ispin,ikpt,ncnt]=np.array([ig1p,ig2p,ig3p]).dot(b)
           # print(igall[ispin,ikpt],ncnt)
            print(time.time()-t1, "igall done")
            wavefile.seek(8*(Recl-4-3*N_Band),1)
            t1=time.time()
            ib0=1
            wavefile.seek(8*Recl*(selectband[0]-ib0),1)
            nb=(selectband[1]-selectband[0]+1)
            coeff_tmp=struct.unpack('f'*2*nb*Recl,wavefile.read(8*nb*Recl))
            coeff_tmp=np.array(coeff_tmp).reshape((nb,Recl,2))
            coeff[ispin,ikpt,selectband[0]-1:selectband[1]]=coeff_tmp[:,:,0]+coeff_tmp[:,:,1]*1j
            norm_coeff=np.linalg.norm(coeff[ispin,ikpt,selectband[0]-1:selectband[1]],axis=1)
            coeff[ispin,ikpt,selectband[0]-1:selectband[1]]=(coeff[ispin,ikpt,selectband[0]-1:selectband[1]].T/norm_coeff).T
            wavefile.seek(8*Recl*(N_Band-selectband[1]),1)
            print(time.time()-t1,"coeff done")
            print("coeff",coeff[0,0,selectband[0]-1])
    wavefile.close()
    print("coeff",coeff[0,0,selectband[0]-1])
    return coeff,igall,eig,occ



def calc_tdm(coeff_,eig_,igall_):
    print('start')
    t1=time.time()
    dE=np.array(np.matrix(eig_)-np.matrix(eig_).T)/(2*13.605826)
    dE[dE==0]=np.inf 
    tdm=np.zeros((4,dE.shape[0],dE.shape[1]),dtype='complex')
    for i in  range(3):
        tdm_=((np.conj(coeff_)*np.array(igall_)[:,i])).dot(coeff_.T)
        tdm_=1j/np.array(dE)*np.array(tdm_)*2.541746*0.529177249
        tdm_=tdm_*np.conj(tdm_)
        tdm[i]=tdm_
    tdm[-1]=np.sum(np.array(tdm),axis=0)
    tdm=np.abs(tdm)
    print('done')
    print(time.time()-t1)
    return tdm

def plot_lamda(lamda,tdm,index,vb,selectband):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax=plt.subplot(111)
    intx=int((selectband[1]-selectband[0])/5)
    if intx <1:
        intx=1
    intx=5
    x=np.arange(selectband[0],selectband[0]+vb+1,intx)
    y=np.arange(selectband[1]+1,selectband[0]+vb+1,intx)
    #x=np.arange(selectband[0],selectband[1]+1,intx)
    #y=np.arange(selectband[0],selectband[1]+1,intx)
    #sc=ax.imshow(tdm[index],origin='lower',cmap="hot_r")
    sns.set_style("dark")
    sns.heatmap(np.abs(lamda[vb+1:,:vb+1]),ax=ax,square=True,cmap='Spectral_r',mask=((np.abs(lamda[vb+1:,:vb+1])>820)|(np.abs(lamda[vb+1:,:vb+1])<350)),vmin=350,vmax=820)
    plt.xticks(list(x-selectband[0]),list(x) ,fontsize=24,rotation='-90')
    plt.yticks(list(y-selectband[0]),list(y),fontsize=24,rotation='0')
    #plt.xticks(list(x-selectband[0]),["VB-%d" % i for i in x[-1]-x[0:-1]]+["VB"] ,fontsize=24,rotation='-90')#-selectaaband[0])
    #plt.yticks(list(y-1-vb-selectband[0]),["CB"]+["CB+%d" % i for i in y[1:]-selectband[0]],fontsize=24,rotation='0')#-selectband[0]) 
    ax.invert_yaxis()
    ax.invert_xaxis()
    #dv=1
    #ax.plot([vb+dv,vb+dv],[vb-dv,tdm.shape[1]-dv],color='b')
    #ax.plot([0-dv,vb+dv],[vb+dv,vb+dv],color='g')
    #ax.plot([vb+dv,tdm.shape[1]-dv],[vb+dv,vb+dv],color='b')
    #ax.plot([vb+dv,vb+dv],[0-dv,vb+dv],color='g')
    #cbar_ax = inset_axes(ax,
    #            width=0.4, # width = 10% of parent_bbox width
    #            height="100%", # height : 50%
    #            loc=6,
    #            bbox_to_anchor=(1.05, 0., 1, 1),
    #            bbox_transform=ax.transAxes,
    #            borderpad=0,
    #        )
    #cb=plt.colorbar(sc,cax=cbar_ax,orientation='vertical')
    #cb.ax.tick_params(labelsize=24)
    print(x,selectband[0])


def plot_tdm(tdm,index,vb,selectband):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax=plt.subplot(111)
    intx=int((vb)/5)
    if intx <1:
        intx=1
    x=np.arange(selectband[0],selectband[0]+vb+1,intx)
    print(x,selectband[0]+vb+1)
    x=x+(selectband[0]+vb-np.max(x))
    y=np.arange(selectband[0]+vb+1,selectband[1]+1,intx)
    sc=ax.imshow(tdm[index,vb+1:,:vb+1],origin='lower',cmap="hot_r")
    plt.xticks(list(x-selectband[0])[:],list(x)[:],fontsize=24,rotation='-90')
    plt.yticks(list(y-1-selectband[0]-vb)[:],list(y)[:],fontsize=24,rotation='0')
    ax.invert_xaxis()
    cbar_ax = inset_axes(ax,
                width=0.4, # width = 10% of parent_bbox width
                height="100%", # height : 50%
                loc=6,
                bbox_to_anchor=(1.05, 0., 1, 1),
                bbox_transform=ax.transAxes,
                borderpad=0,
            )
    cb=plt.colorbar(sc,cax=cbar_ax,orientation='vertical')
    cbar_ax.xaxis.set_label_position('top')
    cbar_ax.set_xlabel(r'$\mathrm{Debye^2}$',fontsize=24)
    cb.ax.tick_params(labelsize=24)
    print(x,selectband[0])

def plot_tdm_byEnergy(Eig,ispin,ikpt,tdm,index,vb,selectband):
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax=plt.subplot(111)
    Eig=Eig[ispin,ikpt,selectband[0]-1:selectband[1]]
    Eintx=(np.max(Eig)-np.min(Eig))/1000
    Emin=np.min(Eig)
    Emax=np.max(Eig)
    EE=np.zeros((1000,1000))
    for i in range(1000):
        for j in range(1000):
            EE[i,j]=tdm[-1,np.sum(Eig<(Emin+i*Eintx))-1,np.sum(Eig<(Emin+j*Eintx))-1]
    x=np.arange(-0.5,999.5,50)*Eintx+Emin
    y=np.arange(-0.5,999.5,50)*Eintx+Emin
#    plt.imshow(EE)
    sns.set_style("whitegrid")
    ax=sns.heatmap(EE,ax=ax,square=True)
    sns.boxplot()
    ax.invert_yaxis()
    plt.xticks(range(0,1000,50),["%2.2f" % i for i in x],fontsize=24,rotation='90')#-selectaaband[0])
    plt.yticks(range(0,1000,50),["%2.2f" % i for i in y],fontsize=24,rotation='0')#-selectband[0])
    #ax.plot(,[vb-dv,tdm.shape[1]-dv],color='b')
    #ax.plot([0-dv,vb+dv],[vb+dv,vb+dv],color='g')
    #ax.plot([vb+dv,tdm.shape[1]-dv],[vb+dv,vb+dv],color='b')
    #ax.plot([vb+dv,vb+dv],[0-dv,vb+dv],color='g')
    #cbar_ax = inset_axes(ax,
    #            width=0.4, # width = 10% of parent_bbox width
    #            height="100%", # height : 50%
    #            loc=6,
    #            bbox_to_anchor=(1.05, 0., 1, 1),
    #            bbox_transform=ax.transAxes,
    #            borderpad=0,
    #        )
    #cb=plt.colorbar(sc,cax=cbar_ax,orientation='vertical')
    #cb.ax.tick_params(labelsize=24)
    print(x,selectband[0])


read=input("read from old data\n")
figsize=(10,5)
#read="T"
if not read is "":
    tdm,iband,ikpt,ispin,vb,eig=np.load(read,allow_pickle=True)
    #tdm,iband,ikpt,ispin,vb,eig=np.load("tdm_0_0_1450_1550.npy",allow_pickle=True)
    tdm[:,:vb+1,:vb+1]=0
    tdm[:,vb+1:,vb+1:]=0
else:
    tmp=(input("Start band,End band,Kpoints,ispin\n")).split()
    print(tmp)
    iband=[int(tmp[0]),int(tmp[1])]
    ikpt=int(tmp[2])
    ispin=int(tmp[3])
    coeff,igall,eig,occ=get_wave(path,iband)
    vb=int(np.sum(occ[ispin,ikpt]>0.001)-iband[0])
    np.save("occ.npy",occ)
    coeff_=coeff[ispin,ikpt]
    eig_=eig[ispin,ikpt]
    igall_=igall[ispin,ikpt]
    tdm=calc_tdm(coeff_[iband[0]-1:iband[1]],eig_[iband[0]-1:iband[1]],igall_)
    np.save("tdm_%d_%d_%d_%d.npy"%(ispin,ikpt,iband[0],iband[1]),(tdm,iband,ikpt,ispin,vb,eig))
    np.savetxt("tdm_"+str(ispin)+'_'+str(iband[0])+'_'+str(iband[1]),tdm[-1])
    #tdm[:,:vb+1,:vb+1]=0
    #tdm[:,vb+1:,vb+1:]=0
eig_=np.array(eig[ispin,ikpt,iband[0]-1:iband[1]]).reshape(iband[1]-iband[0]+1,1)
#lamda=1240/(eig_.T-eig_)
#lamda[tdm[-1]<1]=np.inf
#print(eig_.shape,lamda.shape,(eig_.T-eig_).shape)
#print(iband)
#plot_tdm_byEnergy(eig,ispin,ikpt,tdm,-1,vb,iband)
plt.rc('font',family='Times New Roman')
matplotlib.rcParams['mathtext.fontset'] = 'stix'

fig=plt.figure(figsize=figsize)
plot_tdm(tdm,-1,vb,iband)
plt.tight_layout()
plt.subplots_adjust(right=0.8)

plt.savefig("tdm_%d_%d_%d_%d.png"%(ispin,ikpt,iband[0],iband[1]),dpi=360)
#plt.show()
#fig=plt.figure(figsize=figsize)
#plot_lamda(lamda,tdm,-1,vb,iband)
#plt.tight_layout()
#plt.savefig("lamda_%d_%d_%d_%d.pnb"%(ispin,ikpt,iband[0],iband[1]),dpi=360)
#plt.show()
