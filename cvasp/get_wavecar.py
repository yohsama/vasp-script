import numpy as np
import struct
import scipy
import scipy.integrate
import time
import ase.io
from multiprocessing import Pool
class get_wavecar:
    def get_eig(WAVECARs):
        i=0
        for WAVECAR in WAVECARs:
            wavefile=open(WAVECAR,'rb')
            Recl=int(struct.unpack('d',wavefile.read(8))[0]/8)
            Ispin=int(struct.unpack('d',wavefile.read(8))[0])
            RFLAG=int(struct.unpack('d',wavefile.read(8))[0])
            wavefile.seek(8*(Recl),0)
            N_wave_k_vector=int(struct.unpack('d',wavefile.read(8))[0])
            N_Band=int(struct.unpack('d',wavefile.read(8))[0])
            #coeff=np.zeros((Ispin,N_wave_k_vector,N_Band,Recl),dtype='double')
            #eig=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
            eig_=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
            occ_=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
            wave_k_vector_=np.zeros((N_wave_k_vector,3),dtype='double')
            for ispin in range(Ispin):
                for ikpt in range(N_wave_k_vector):
                    wavefile.seek(8*Recl*(2+((N_wave_k_vector*ispin+ikpt)*(1+N_Band))),0)
                    wavefile.read(8)
                    for ik in range(3):
                        wave_k_vector_[ikpt,ik]=struct.unpack('d',wavefile.read(8))[0]
                    #t1=time.time()
                    eig_tmp=np.array(struct.unpack('d'*3*N_Band,wavefile.read(8*3*N_Band)));
                    eig_[ispin,ikpt]=np.sign(eig_tmp[0::3])*np.abs(eig_tmp[0::3]+eig_tmp[1::3]*1j)
                    occ_[ispin,ikpt]=eig_tmp[2::3]
                    wavefile.seek(8*Recl*(2+1+((N_wave_k_vector)*(ispin)+ikpt)*(1+N_Band)),0)
                   #t1=time.time()
            wavefile.close()
            if i==0:
                eig=eig_
                wave_k_vector=wave_k_vector_
                occ=occ_
                i=1
            else:
                eig=np.concatenate((eig,eig_),axis=1)
                wave_k_vector=np.concatenate((wave_k_vector,wave_k_vector_),axis=0)
                occ=np.concatenate((occ,occ_),axis=1)
        return eig,occ
    
    def get_igall(WAVECARs,CONTCAR):
        #t1=time.time()
        c=0.26246582250210965422
        try:
            SYS=ase.io.read(CONTCAR,format='vasp')
        except:
            print("can't find CONTCAR")
            return
        i=0
        for WAVECAR in WAVECARs:
            wavefile=open(WAVECAR,'rb')
            Recl=int(struct.unpack('d',wavefile.read(8))[0]/8)
            Ispin=int(struct.unpack('d',wavefile.read(8))[0])
            RFLAG=int(struct.unpack('d',wavefile.read(8))[0])
            wavefile.seek(8*(Recl),0)
            N_wave_k_vector=int(struct.unpack('d',wavefile.read(8))[0])
            N_Band=int(struct.unpack('d',wavefile.read(8))[0])
            Encut=(struct.unpack('d',wavefile.read(8))[0])
            cell=np.zeros((3,3))
            for ci in range(3):
                for cj in range(3):
                    cell[ci,cj]=struct.unpack('d',wavefile.read(8))[0]
            b=np.zeros((3,3))
            Vcell=np.dot(np.cross(cell[0,:],cell[1,:]),cell[2,:])
            b[0,:]=2*np.pi*np.cross(cell[1,:],cell[2,:])/Vcell
            b[1,:]=2*np.pi*np.cross(cell[2,:],cell[0,:])/Vcell
            b[2,:]=2*np.pi*np.cross(cell[0,:],cell[1,:])/Vcell
            bmag=np.sqrt(np.sum(b**2,1))
            #phi12=np.arccos(np.dot(b[0,:],b[1,:])/bmag[0]*bmag[1])
            phi12=np.arccos(np.dot(b[0,:],b[1,:])/bmag[0]/bmag[1])
            #print(np.dot(b[0,:],b[1,:]))
            vtmp=np.cross(b[0,:],b[1,:])
            vmag=np.linalg.norm(np.cross(b[0,:],b[1,:]))
            sinphi123=np.dot(b[2,:],vtmp)/(vmag*bmag[2])
            nb1maxA=np.round(np.sqrt(Encut*c)/(bmag[0]*abs(np.sin(phi12))))+1
            nb2maxA=np.round(np.sqrt(Encut*c)/(bmag[1]*abs(np.sin(phi12))))+1
            nb3maxA=np.round(np.sqrt(Encut*c)/(bmag[2]*abs(sinphi123)))+1
            npmaxA=np.round(4*np.pi*nb1maxA*nb2maxA*nb3maxA/3)
            #phi13=np.arccos(np.dot(b[0,:],b[2,:])/bmag[0]*bmag[2])
            phi13=np.arccos(np.dot(b[0,:],b[2,:])/bmag[0]/bmag[2])
            vtmp=np.cross(b[0,:],b[2,:])
            vmag=np.linalg.norm(vtmp)
            sinphi123=np.dot(b[1,:],vtmp)/(vmag*bmag[1])
            nb1maxB=np.round(np.sqrt(Encut*c)/bmag[0]*abs(np.sin(phi13))+1)
            nb2maxB=np.round(np.sqrt(Encut*c)/bmag[1]*abs(sinphi123)+1)
            nb3maxB=np.round(np.sqrt(Encut*c)/bmag[2]*abs(np.sin(phi13))+1)
            npmaxB=np.round(4*np.pi*nb1maxB*nb2maxB*nb3maxB/3)
            #phi23=np.arccos(np.dot(b[1,:],b[2,:])/bmag[1]*bmag[2])
            phi23=np.arccos(np.dot(b[1,:],b[2,:])/bmag[1]/bmag[2])
            vtmp=np.cross(b[1,:],b[2,:])
            vmag=np.linalg.norm(vtmp)
            sinphi123=np.dot(b[0,:],vtmp)/(vmag*bmag[1])
            nb1maxC=np.round(np.sqrt(Encut*c)/bmag[0]*abs(sinphi123)+1)
            nb2maxC=np.round(np.sqrt(Encut*c)/bmag[1]*abs(np.sin(phi23))+1)
            nb3maxC=np.round(np.sqrt(Encut*c)/bmag[2]*abs(np.sin(phi23))+1)
            npmaxC=np.round(4*np.pi*nb1maxC*nb2maxC*nb3maxC/3)
            nb1max=int(np.max([nb1maxA,nb1maxB,nb1maxC]))
            nb2max=int(np.max([nb2maxA,nb2maxB,nb2maxC]))
            nb3max=int(np.max([nb3maxA,nb3maxB,nb3maxC]))
            npmax=int(np.max([npmaxA,npmaxB,npmaxC]))
            sumkg=np.zeros(3)
            wave_k_vector_=np.zeros((N_wave_k_vector,3),dtype='double')
            igall_=np.zeros((Ispin,N_wave_k_vector,Recl,3),dtype='double')
            nb1=list(range(0,nb1max+1))+list(range(-nb1max,0))
            nb2=list(range(0,nb2max+1))+list(range(-nb2max,0))
            nb3=list(range(0,nb3max+1))+list(range(-nb3max,0))
            igtmp=np.array(np.meshgrid(nb2,nb3,nb1),dtype='int').reshape(3,(2*nb3max+1)*(2*nb2max+1)*(2*nb1max+1)).T
            igtmp=igtmp[:,np.array([2,0,1])]
            for ispin in range(Ispin):
                for ikpt in range(N_wave_k_vector):
                    wavefile.seek(8*Recl*(2+(N_wave_k_vector*ispin+ikpt)*(N_Band+1)),0)
                    N_plane=int(struct.unpack('d',wavefile.read(8))[0])
                    for ik in range(3):
                        wave_k_vector_[ikpt,ik]=struct.unpack('d',wavefile.read(8))[0]
                    #igtmp=igtmp[:,np.array([2,0,1])]
#np.savetxt('igtmp',igtmp)
                    sumkg=np.matrix(igtmp+wave_k_vector_[ikpt])*np.matrix(b)
                    gtot=np.linalg.norm(sumkg,axis=1)
                    etot=(gtot**2)/c
                    igtmp_=igtmp[etot<=Encut]
                    #igtmp_=np.matrix(igtmp_)*np.matrix(b)
                    igall_[ispin,ikpt,:igtmp_.shape[0],:]=igtmp_
            if i==0:
                igall=igall_
                i=1
            else:
                if igall.shape[2]<igall_.shape[2]:
                    tmp1=np.zeros((igall.shape[0],igall.shape[1],igall_.shape[2],igall.shape[3]))
                    tmp1[:,:,:igall.shape[2]]=igall
                    igall=np.concatenate((tmp1,igall_),axis=1)
                else:
                    tmp2=np.zeros((igall_.shape[0],igall_.shape[1],np.max((igall.shape[2],igall_.shape[2])),igall_.shape[3]))
                    tmp2[:,:,:igall_.shape[2]]=igall_
                    igall=np.concatenate((igall,tmp2),axis=1)
        #print(igall)
        return igall,b
    
    
    
    def get_coeff(WAVECARs,iband):
        i=0
        nb=len(iband)
        Nkpt=0
        MRecl=0  
        for WAVECAR in WAVECARs:
            wavefile=open(WAVECAR,'rb')
            Recl=int(struct.unpack('d',wavefile.read(8))[0]/8)
            Ispin=int(struct.unpack('d',wavefile.read(8))[0])
            RFLAG=int(struct.unpack('d',wavefile.read(8))[0])
            wavefile.seek(8*(Recl),0)
            N_wave_k_vector=int(struct.unpack('d',wavefile.read(8))[0])
            N_Band=int(struct.unpack('d',wavefile.read(8))[0])
            c=0.26246582250210965422
            Nkpt=Nkpt+N_wave_k_vector
            MRecl=np.max((Recl,MRecl))
        coeff=np.zeros((Ispin,Nkpt,N_Band,MRecl),dtype='complex')
        iNkpt=0
        for WAVECAR in WAVECARs:
            t1=time.time()
            wavefile=open(WAVECAR,'rb')
            Recl=int(struct.unpack('d',wavefile.read(8))[0]/8)
            Ispin=int(struct.unpack('d',wavefile.read(8))[0])
            RFLAG=int(struct.unpack('d',wavefile.read(8))[0])
            wavefile.seek(8*(Recl),0)
            N_wave_k_vector=int(struct.unpack('d',wavefile.read(8))[0])
            N_Band=int(struct.unpack('d',wavefile.read(8))[0])
            c=0.26246582250210965422
            #coeff=np.zeros((Ispin,N_wave_k_vector,N_Band,Recl),dtype='double')
            #eig=np.zeros((Ispin,N_wave_k_vector,N_Band),dtype='double')
            for ispin in range(Ispin):
                for ikpt in range (N_wave_k_vector):
                    N=3+((N_wave_k_vector)*(ispin)+ikpt)*(1+N_Band)
                    wavefile.seek(8*Recl*N,0)
                    wavefile.seek(8*Recl*(iband[0]-1),1)
                    coeff_tmp=struct.unpack('f'*2*nb*Recl,wavefile.read(8*nb*Recl))
                    coeff_tmp=np.array(coeff_tmp).reshape((nb,Recl,2))
                    coeff_jtmp=coeff_tmp[:,:,0]+coeff_tmp[:,:,1]*1j
                    coeff[ispin,iNkpt+ikpt,iband-1,:Recl]=coeff_jtmp
            iNkpt=iNkpt+N_wave_k_vector
                    #print(time.time()-t1,"coeff done")
            wavefile.close()
        t1=time.time()
        bcoeff=coeff[:,:,iband-1]
        bcoeff_norm=np.linalg.norm(bcoeff,axis=3).repeat(bcoeff.shape[3]).reshape(bcoeff.shape)
        coeff[:,:,iband-1]=bcoeff/bcoeff_norm
        return coeff
    
