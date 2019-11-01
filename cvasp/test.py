#!/usr/bin/python3
from plot_band import *
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from glob import glob
ISPIN=int(input("ISPIN\n"))
plottype=int(input("fat_band/1,colormap_band/2\n"))
type1_Ele=[]
type2_Ele_A=[]
type2_Ele_B=""
if plottype%10==1:
    tmp=(input("Element,split by \",\" color symbol\n"))
    while not tmp == "":
        type1=tmp.split()
        print(type1)
        Element=type1[0].split(",")
        #if len(type1)<4:
        #    type1.append("")
        type1_Ele.append([Element,type1[1]])
        tmp=(input("Element,split by \",\" color symbol; Enter to end\n"))


if plottype==2:
    tmp=(input("Element,split by \",\"\n"))
    type2_Ele_A=tmp.split(",")
    tmp=input("Part B ? ENTER for EXIT, any words for CONTINUE\n")
    if not tmp =="":
        type2_Ele_B=tmp.split(",")
    else:
        type2_Ele_B=False

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
    project=pd.read_csv('./project.csv')
    print("read successful")
    kpoints=plot_band.from_kpoints(path='.',read=False,save=True)
except Exception as e:
    print("try to read but fail\n",e)
    project=plot_band.from_procar(ISPIN=ISPIN,read=True)
    kpoints=plot_band.from_kpoints(path='.',KPOINTS='KPOINTS',read=False,save=True)
    #project=project.merge(kpoints[['kpath','kpt','kpt_label']])
    #project=plot_band.add_info_symbol(project=project,CONTCAR_file='./CONTCAR')
    project['eig']=project['eig']-project[project['occ']>0.001]['eig'].max()
    project.to_csv('./project.csv')

project=plot_band.add_info_symbol(project=project,CONTCAR_file='./CONTCAR')
for ispin in range(ISPIN):
    fig=plt.figure(figsize=figsize)
    #if plot_type!=2:
    # set special K-points
    vbm=project[project['occ']>0.001]['eig'].max()
    cbm=project[project['occ']<=0.001]['eig'].min()
    dElim=(Elim[1]-Elim[0])
    if False: #(cbm-vbm)>(dElim/10):
        try:
            from brokenaxes import brokenaxes
 #           print(vbm+dElim/20,cbm-dElim/20)
            ax = brokenaxes(ylims=((Elim[0],vbm+dElim/20),(cbm-dElim/20,Elim[1])), subplot_spec=())
            print('yes')
        except Exception as e:
            print(e)
            ax=plt.subplot()
            ax.set_ylim((Elim[0],Elim[1]))
    else:
        ax=plt.subplot()
        ax.set_ylim((Elim[0],Elim[1]))
    #ax.set_xlim([kpoints['kpath'].min(),kpoints['kpath'].max()])
    #ax.set_yticklabels(ax.get_yticks(),fontsize=20)

    ax,ax_cbar,divider=plot_band.plot_pband(ax=ax,plot_type=plottype,kpoints=kpoints,project=project,ispin=ispin,Elim=Elim,type1_Ele=type1_Ele,type2_Ele_A=type2_Ele_A,type2_Ele_B=type2_Ele_B)
    ax=plot_band.plot_sp_kline(kpoints,ax=ax)

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
        ax_tdm=plot_band.plot_sp_kline(kpoints,ax=ax_tdm)
        tdm_k=np.array([])
        ymax=0
        for iselectband in selectband:
            eig,occ=get_wavecar.get_eig(WAVECARs)
            igall,b=get_wavecar.get_igall(WAVECARs,CONTCAR)
            iband=np.arange(iselectband[0],iselectband[1]+1)
            coeff=get_wavecar.get_coeff(WAVECARs,iband)
            tdm_k,tdm=cal_tdm_byband(ispin,iselectband,coeff,igall,eig,occ,b)
            np.savetxt('tdm_cx.dat',tdm)
            print(tdm.shape)
            ax_tdm=plot_tdm_band(kpoints,tdm[:,-1],0,label=iselectband,ax_tdm=ax_tdm)
            ymax=np.max((ymax,np.max(tdm_k)))
        #ax_tdm.set_yticklabels(ax.get_yticks())#,fontsize=24)
        ax_tdm.set_ylim((0,ymax*1.1))
        plt.legend( markerscale=0.5)#,fontsize=24)
    plt.tight_layout()
    plt.savefig('band_%d.png' % ispin,dpi=300)
plt.show()
