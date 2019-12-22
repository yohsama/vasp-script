#!/usr/bin/python3
from plot_band import *
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from glob import glob
from set_group import set_group
ISPIN=int(input("ISPIN\n"))
plottype=int(input("fat_band/1,colormap_band/2\n"))
type1_Ele=[]
type2_Ele_A=[]
type2_Ele_B=""
GROUP=[]
GROUP_INFO=[]
if plottype%10==1:
    tmp=(input("Element,split by ',' color tag \n"))
    while not tmp == "":
        type1=tmp.split()
        #print(type1)
        GROUP.append(type1[0])
        GROUP_INFO.append(type1[1:])
        #if len(type1)<4:
        #    type1.append("")
        tmp=(input("Element,split by \",\" color symbol; Enter to end\n"))


if plottype==2:
    tmp=(input("Element,split by \",\"\n"))
    GROUP.append(tmp.split()[0])
    GROUP_INFO.append(tmp.split()[1])
    tmp=input("Part B ? ENTER for EXIT, any words for CONTINUE\n")
    if not tmp =="":
        type1=tmp.split()
        GROUP.append(type1[0])
        GROUP_INFO.append(type1[1])
    else:
        type2_Ele_B=False
        GROUP_INFO.append("")

tmp=input("Energy ymin ymax\n").split()
print(tmp)
Elim=(float(tmp[0]),float(tmp[1]))

tmp=input("tdm band \n")
selectband=[]
tdm=False
while not tmp =="":
    tdm=True;
    tmp=tmp.split()
    selectband.append([int(tmp[0]),int(tmp[1])])
    tmp=input("tdm band \n")
    print(tmp)
else:
    pass

tmp=input("figsize \n").split()
if not tmp =="":
    figsize=(float(tmp[0]),float(tmp[1]))
else:
    figsize=(6,8)

tmp=input("fontsize \n").split()
if not tmp =="":
    fontsize=float(tmp[0])
else:
    fontsize='24'

plt.rc('font',family='Times New Roman',size=fontsize)
matplotlib.rcParams['mathtext.fontset'] = 'stix'

#if True:
try :
    project=np.load('./project.npy',allow_pickle=True)
    print("read successful")
    kpath=plot_band.from_kpoints(path='.',read=False,save=True)
except Exception as e:
    print("try to read but fail\n",e)
    #PRO,EIG,OCC,TTOT,KPOINTS
    project=plot_band.from_procar(ISPIN=ISPIN,read=True)
    #KPATH,KLABEL,KPOINTS=
    kpath=plot_band.from_kpoints(path='.',KPOINTS='KPOINTS',read=False,save=True) 
    #project=project.merge(kpoints[['kpath','kpt','kpt_label']])
    #project=plot_band.add_info_symbol(project=project,CONTCAR_file='./CONTCAR')
    #project['eig']=project['eig']-project[project['occ']>0.001]['eig'].max()
    np.save('./project.npy',project)

PRO,EIG,OCC,TTOT,KPOINTS=project
KPATH,KLABEL,KPOINTS=kpath
#print(KLABEL)
set_fermi=input('set fermi to zero?')
if 'T' in set_fermi:
    FERMI=plot_band.get_fermi(EIG,OCC)
    #project['eig']=project['eig']-project[project['occ']>0.001]['eig'].max()
elif not set_fermi is "":
    FERMI=np.float(set_fermi)
else:
    FERMI=0
    #project['eig']=project['eig']-float(set_fermi)
SYMBOL=plot_band.get_symbol(CONTCAR_file='./CONTCAR')
for ispin in range(EIG.shape[0]):
    fig=plt.figure(figsize=figsize)
    #if plot_type!=2:
    # set special K-points
    ax=plt.subplot()
    ax.set_ylim((Elim[0],Elim[1]))
    if kpath[0].min()==kpath[0].max():
        ax.set_xlim([-0.5,0.5])
        ax.set_xticklabels([])
    else:
        ax.set_xlim([kpath[0].min(),kpath[0].max()])

    #ax.set_xlim([kpoints['kpath'].min(),kpoints['kpath'].max()])
    #ax.set_yticklabels(ax.get_yticks(),fontsize=20)
    PRO_GROUP=set_group(PRO,SYMBOL,group=GROUP)
    ax,ax_cbar,divider=plot_band.plot_pband(
                    plot_type=plottype,
                    kpath=kpath,
                    eig=EIG[ispin],
                    pro_group=PRO_GROUP[:,ispin],
                    group_info=GROUP_INFO,
                    ax=ax,
                    ispin=ispin,
                    fermi=FERMI,
                    )

    ax=plot_band.plot_sp_kline(kpath,ax=ax)

    if tdm:
        from plot_tdm_band import *
        CONTCAR=path+'/CONTCAR'
        WAVECARs=glob(path+'/WAVECAR_*')
        if WAVECARs ==[] :
            WAVECARs=[path+'/WAVECAR']
        WAVECARs.sort()
        ax.set_xticks=([])
        ax.set_xticklabels([])
        ax_tdm= divider.append_axes('bottom', size="15%", pad=0.1)
        ax_tdm=plot_band.plot_sp_kline(kpath,ax=ax_tdm)
        tdm_k=np.array([])
        ymax=0
        for iselectband in selectband:
            eig,occ=get_wavecar.get_eig(WAVECARs)
            igall,b=get_wavecar.get_igall(WAVECARs)
            iband=np.arange(iselectband[0],iselectband[1]+1)
            coeff=get_wavecar.get_coeff(WAVECARs,iband)
            tdm_k,tdm=cal_tdm_byband(ispin,iselectband,coeff,igall,eig,occ,b)
            np.savetxt('tdm_cx.dat',tdm)
            #print(tdm.shape)
            ax_tdm=plot_tdm_band(kpath,tdm[:,-1],0,label=iselectband,ax_tdm=ax_tdm)
            ymax=np.max((ymax,np.max(tdm_k)))
        #ax_tdm.set_yticklabels(ax.get_yticks())#,fontsize=24)
        ax_tdm.set_ylim((0,ymax*1.1))
        plt.legend( markerscale=0.5)#,fontsize=24)
    plt.tight_layout()
    plt.savefig('band_%d.png' % ispin,dpi=300)
plt.show()
