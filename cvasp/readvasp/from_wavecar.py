import struct
import time
from multiprocessing import Pool

import ase.io
import numpy as np
import scipy
import scipy.integrate
from constant import *



class get_wavecar(object):
    def __init__(self,WAVECAR='./WAVECAR',L_gam=False,L_gamma_half='x',L_soc=False):
        self.L_gam=L_gam
        self.L_soc=L_soc
        self._L_gamma_half_=L_gamma_half
        self.wavecar=WAVECAR
        self._file_=open(self.wavecar,'rb')
        self._Recl_,self.I_spin,self._Rflag_=np.array(np.fromfile(self._file_,dtype=np.double,count=3),dtype=int)
        self._Len_coeff_=self._Recl_//8
        self._prec_coeff_=self.get_coeff_prec()
        self._file_.seek((self._Recl_),0)
        self.N_kpt,self.N_Band,self.Encut=np.array(np.fromfile(self._file_,dtype=np.double,count=3),dtype=int)
        self.cell=self.get_cell()
        self.bcell=self.rcell_to_bcell(self.cell)
        self.kpoint=np.zeros((self.N_kpt,3),dtype='double')
        self.eig,self.occ=self.get_eig_and_occ()
        self.I_gvector,self.gvector=self.get_I_gvector(init=True)
        if self._len_gvector_[0] != self._N_planes_[0]:
            if self._len_gvector_[0] == self._N_planes_[0]*2-1 :
                print("Plane number not equal to the igtmp, now rerun this model as Gamma version calculaion.")
                self.L_gam=True
                self.I_gvector,self.gvector=self.get_I_gvector()
            if self._len_gvector_[0] == self._N_planes_[0]/2 :
                print("Plane number not equal to the igtmp, now rerun this model as SOC version calculaion.")
                self.L_soc=True
                self.I_gvector,self.gvector=self.get_I_gvector()
        else:
            self.I_gvector,self.gvector=self.get_I_gvector()
        self._shapes_=np.array((self.I_spin,self.N_kpt,self.N_Band,self._Len_coeff_))
        self.coeff=np.zeros(self._shapes_,self._prec_coeff_)
        self.coeff.reshape(self._shapes_)
        self.tdm=get_tdm(self)

    def get_coeff_prec(self):
        '''
            Rtag = 45200: single precision complex
            Rtag = 45210: double precision complex
        '''
        if self._Rflag_ == 45200:
            return np.complex64
        elif self._Rflag_ == 45210:
            return np.complex128
        else:
            raise ValueError("Invalid TAG values: {}".format(self._Rflag_))


    def get_cell(self):
        cell=np.fromfile(self._file_,np.double,count=9).reshape((3,3))
        return cell
    
    def get_eig_and_occ(self):
        eig=np.zeros((self.I_spin,self.N_kpt,self.N_Band),dtype='double')
        occ=np.zeros((self.I_spin,self.N_kpt,self.N_Band),dtype='double')
        self._N_planes_=np.zeros((self.N_kpt),dtype=int)
        for ispin in range(self.I_spin):
            for ikpt in range(self.N_kpt):
                self._file_.seek(self._Recl_*(2+((self.N_kpt*ispin+ikpt)*(1+self.N_Band))),0)
                self._N_planes_[ikpt]=np.fromfile(self._file_,np.double,count=1)
                self.kpoint[ikpt]=np.fromfile(self._file_,np.double,count=3)
                eig_tmp=np.fromfile(self._file_,np.double,count=(self.N_Band*3)).reshape(-1,3)
                """
                eig_tmp
                            eigenvalue_real            eigenvalue_image              occupid               
                                xxxxx                       xxxxx                      0-2

                """
                eig[ispin,ikpt]=np.sign(eig_tmp[:,0])*np.abs(eig_tmp[:,0]+eig_tmp[:,1]*1j)
                occ[ispin,ikpt]=eig_tmp[:,2]
        return eig,occ

    def rcell_to_bcell(self,cell):
        bcell=np.linalg.inv(cell).T*2*np.pi
        return bcell


    def get_fftgrid(self):
        bmag=np.linalg.norm(self.bcell,axis=1)
        phi12=np.arccos(np.dot(self.bcell[0],self.bcell[1])/bmag[0]/bmag[1])
        vtmp=np.cross(self.bcell[0],self.bcell[1])
        vmag=np.linalg.norm(vtmp)
        sinphi123=np.dot(self.bcell[2],vtmp)/(vmag*bmag[2])
        nb1maxA=np.sqrt(self.Encut*c)//(bmag[0]*abs(np.sin(phi12)))
        nb2maxA=np.sqrt(self.Encut*c)//(bmag[1]*abs(np.sin(phi12)))
        nb3maxA=np.sqrt(self.Encut*c)//(bmag[2]*abs(sinphi123))
        phi13=np.arccos(np.dot(self.bcell[0],self.bcell[2])/bmag[0]/bmag[2])
        vtmp=np.cross(self.bcell[0],self.bcell[2])
        vmag=np.linalg.norm(vtmp)
        sinphi123=np.dot(self.bcell[1],vtmp)/(vmag*bmag[1])
        nb1maxB=np.sqrt(self.Encut*c)//bmag[0]*abs(np.sin(phi13))
        nb2maxB=np.sqrt(self.Encut*c)//bmag[1]*abs(sinphi123)
        nb3maxB=np.sqrt(self.Encut*c)//bmag[2]*abs(np.sin(phi13))
        phi23=np.arccos(np.dot(self.bcell[1],self.bcell[2])/bmag[1]/bmag[2])
        vtmp=np.cross(self.bcell[1],self.bcell[2])
        vmag=np.linalg.norm(vtmp)
        sinphi123=np.dot(self.bcell[0],vtmp)/(vmag*bmag[1])
        nb1maxC=np.sqrt(self.Encut*c)//bmag[0]*abs(sinphi123)
        nb2maxC=np.sqrt(self.Encut*c)//bmag[1]*abs(np.sin(phi23))
        nb3maxC=np.sqrt(self.Encut*c)//bmag[2]*abs(np.sin(phi23))
        nb1max=int(np.max([nb1maxA,nb1maxB,nb1maxC]))
        nb2max=int(np.max([nb2maxA,nb2maxB,nb2maxC]))
        nb3max=int(np.max([nb3maxA,nb3maxB,nb3maxC]))
        nb1=list(range(0,nb1max+1))+list(range(-nb1max,0))
        nb2=list(range(0,nb2max+1))+list(range(-nb2max,0))
        nb3=list(range(0,nb3max+1))+list(range(-nb3max,0))
        return nb1,nb2,nb3

    def get_I_gvector(self,init=False):
        #t1=time.time()
        nb1,nb2,nb3=self.get_fftgrid()
        igtmp=np.array(np.meshgrid(nb2,nb3,nb1),dtype='int').reshape(3,-1).T
        igtmp=igtmp[:,np.array([2,0,1])]
        sumkg=np.zeros(3)
        I_gvector=np.zeros((self.N_kpt,self._Len_coeff_,3),dtype=np.int)
        gvector=np.zeros((self.N_kpt,self._Len_coeff_,3),dtype=np.float)
        self._len_gvector_=np.zeros((self.N_kpt),dtype=int)
        if self.L_gam:
            if self._L_gamma_half_ == 'x':
                print('x')
                igtmp = igtmp[
                    (igtmp[:,0]>0)  |   ((igtmp[:,0]==0) & (igtmp[:,1]>0))  |   ((igtmp[:,0]==0) & (igtmp[:,1]==0) & (igtmp[:,2]>=0))
                    ]
            else:
                print('z')
                igtmp = igtmp[
                    (igtmp[:,2]>0)  |   ((igtmp[:,2]==0) & (igtmp[:,1]>0))  |   ((igtmp[:,2]==0) & (igtmp[:,1]==0) & (igtmp[:,0]>=0))
                    ]
        for ikpt in range(self.N_kpt):
            sumkg=np.matrix(igtmp+self.kpoint[ikpt])*np.matrix(self.bcell)
            gtot=np.linalg.norm(sumkg,axis=1)
            etot=(gtot**2)/c
            igtmp=igtmp[etot<=self.Encut]
            self._len_gvector_[ikpt]=igtmp.shape[0]
            if not init:
                I_gvector[ikpt,:igtmp.shape[0]]=igtmp
                gvector[ikpt,:igtmp.shape[0]]=igtmp.dot(self.bcell)
        return I_gvector,gvector

    def get_single_coeff(self,ispin=1,ikpt=1,iband=1,norm=False):
        """
        read the single point wavefunction coefficient, and normlized
        default:
        ispin=1,ikpt=1,iband=1,norm=False
        """
        N=3+(self.N_kpt*(ispin-1)+(ikpt-1))*(1+self.N_Band)+(iband-1)
        self._file_.seek(self._Recl_*N,0)
        coeff_tmp=np.fromfile(self._file_,dtype=self._prec_coeff_,count=self._Len_coeff_)
        if norm:
            coeff_tmp /= np.linalg.norm(coeff_tmp)
        return coeff_tmp
    
    def get_all_coeff(self,norm=False):
        """
        read the single point wavefunction coefficient, and normlized
        default:
        norm=False
        """
        for ispin in range(self.I_spin):
            for ikpt in range(self.N_kpt):
                N=3+((self.N_kpt)*(ispin)+ikpt)*(1+self.N_Band)
                self._file_.seek(self._Recl_*N,0)
                coeff_tmp=np.fromfile(self._file_,dtype=self._prec_coeff_,count=self.N_Band*self._Len_coeff_).reshape((self.N_Band,self._Len_coeff_))
                self.coeff[ispin,ikpt]=coeff_tmp
        if norm:
            self.coeff = self.coeff/np.linalg.norm(self.coeff,axis=3)[:,:,:,np.newaxis]
        return self.coeff

    def get_tdm(self):
        self.tdm=get_tdm(self)
        return self.tdm

class get_tdm(object):
    def __init__(self,wfc):
        self.wfc=wfc
        pass

    def calc_single_tdm(self,ikpt,ispin1,iband1,ispin2,iband2,norm=True):
        """
        this module will calculate the Transition dipole moment 
        between two KS-band
        according to the formula in:
        https://en.wikipedia.org/wiki/Transition_dipole_moment
                                            i⋅h
        <psi_a|r|psi_b> = ( i * hbar) / ( m ( Eb - Ea )) * <psi_a|p|psi_b>
                        = ( hbar )^2 / ( m ( Eb - Ea )) * SUM_i( Cai * Cbi *Gi )
        where:
            a, b for the state of KS-band
            C for the coefficient of KS-band
            G for the gvector of specific kpoint

        input as
        ikpt=1, 
        ispin1=1, rband1 = 12
        ispin1=2, rband2 = 33
        norm    = True      this option determines whether to normalize the coefficients 

        returns Transition dipole moment and overlap

        """
        wfc=self.wfc
        dE=wfc.eig[ispin1-1,ikpt-1,iband1-1]-wfc.eig[ispin2-1,ikpt-1,iband2-1]
        coeff1=wfc.get_single_coeff(ispin1,ikpt,iband1,norm=norm)
        igall=wfc.gvector[ikpt-1]
        coeff2=wfc.get_single_coeff(ispin2,ikpt,iband2,norm=norm)
        tdm=np.zeros(4,dtype=np.complex128)
        tmp=coeff1.conj()*coeff2
        overlap=tmp.sum()
        tmp_project_g=np.sum(tmp[:,np.newaxis]*igall,axis=0)
        if wfc.L_gam:
            tmp2=coeff1*coeff2.conj()
            tmp_project_g_r =   np.sum(tmp2[:,np.newaxis]*igall,axis=0)
            tmp_project_g   =  (tmp_project_g-tmp_project_g_r)/2
        tdm_=1j/dE*(2*hartree2eV)*tmp_project_g*au2debye*au2anstrom
        tdm[:3]=tdm_
        tdm[-1]=np.sum(tdm_*tdm_.conj(),axis=0)
        return tdm , overlap

    def calc_tdm_in_kpoint(self,ikpt,ispin1,rband1,ispin2,rband2,norm=True):
        """
        his module will calculate the Transition dipole moment 
        between multiple KS-band in specific kpoint
        according to the formula in:
        https://en.wikipedia.org/wiki/Transition_dipole_moment
                                            i⋅h
        <psi_a|r|psi_b> = ( i * hbar) / ( m ( Eb - Ea )) * <psi_a|p|psi_b>
                        = ( hbar )^2 / ( m ( Eb - Ea )) * SUM_i( Cai * Cbi *Gi )
        where:
            a, b for the state of KS-band
            C for the coefficient of KS-band
            G for the gvector of specific kpoint

        input as:

        ikpt=1, 
        ispin1=1, rband1 = np.array([1,2,3]
        ispin1=2, rband2 = np.array([1,2,3]
        norm    = True      this option determines whether to normalize the coefficients 

        returns TDM and overlap as :

                band1_1 band1_2 band1_3 band1_4
        band2_1 xxxx    xxxx    xxxx    xxxx
        band2_2 xxxx    xxxx    xxxx    xxxx    
        band2_3 ...     ... 
        """
        wfc=self.wfc
        band_index1=rband1.shape[0]
        band_index2=rband2.shape[0] 
        dE=(wfc.eig[ispin1-1,ikpt-1,rband1-1]-(wfc.eig[ispin2-1,ikpt-1,rband2-1])[:,np.newaxis])
        igall=wfc.gvector[ikpt-1]
        coeff1=np.zeros((band_index1,igall.shape[0]),np.complex)
        coeff2=np.zeros((band_index2,igall.shape[0]),np.complex)
        for i,iband in enumerate(rband1):
            coeff1[i]=wfc.get_single_coeff(ispin1,ikpt,iband,norm=norm)
        for i,iband in enumerate(rband2):
            coeff2[i]=wfc.get_single_coeff(ispin2,ikpt,iband,norm=norm)
        tdm=np.zeros((band_index1,band_index2,4),dtype=np.complex128)
        tmp=coeff1[np.newaxis].conj()*coeff2[:,np.newaxis]
        overlap=tmp.sum(2)
        tmp_project_g=np.sum(tmp[:,:,:,np.newaxis]*igall[np.newaxis,np.newaxis],axis=2)
        if wfc.L_gam:
            tmp2=coeff1[np.newaxis]*coeff2[:,np.newaxis].conj()
            tmp_project_g_r =   np.sum(tmp2[:,:,:,np.newaxis]*igall[np.newaxis,np.newaxis],axis=2)
            tmp_project_g   =  (tmp_project_g-tmp_project_g_r)/2
        tdm_=1j/dE[:,:,np.newaxis]*(2*hartree2eV)*tmp_project_g*au2debye*au2anstrom
        tdm[:,:,:3]=tdm_
        tdm[:,:,-1]=np.sum(tdm_*tdm_.conj(),axis=-1)
        #tdm=np.abs(tdm)
        return tdm,overlap

    def calc_tdm_in_band(self,kpath,ispin1,iband1,ispin2,iband2,norm=True):
        """
        his module will calculate the Transition dipole moment 
        between multiple KS-band in specific kpoint
        according to the formula in:
        https://en.wikipedia.org/wiki/Transition_dipole_moment
                                            i⋅h
        <psi_a|r|psi_b> = ( i * hbar) / ( m ( Eb - Ea )) * <psi_a|p|psi_b>
                        = ( hbar )^2 / ( m ( Eb - Ea )) * SUM_i( Cai * Cbi *Gi )
        where:
            a, b for the state of KS-band
            C for the coefficient of KS-band
            G for the gvector of specific kpoint
        
        input as:

        kpath=np.array( [ kpt1, kpt2, kpt3, ...] ), 
        ispin1  =1, rband1 = 23
        ispin1  =2, rband2 = 11
        norm    = True      this option determines whether to normalize the coefficients 
        returns TDM and overlap as :

                tdm
        kpoint1 xxxx   
        kpoint2 xxxx   
        kpoint3 ...     

        """
        tdm=np.zeros((kpath.size,4),dtype=np.complex)
        overlap=np.zeros((kpath.size),dtype=np.complex)
        for ikpt in kpath:
            #print("handling kpoint %d, Total is %d" % (ikpt,len(kpath)))
            tdm[ikpt-1],overlap[ikpt-1]=self.calc_single_tdm(ikpt,ispin1,iband1,ispin2,iband2,norm=norm)
        return tdm  

if __name__ == '__main__':
 
    pass
